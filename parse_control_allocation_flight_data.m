%% Parse PX4 control-allocation replay data
% 这个脚本只做一件事：把 .ulg 里的 PX4 控制分配输入/输出整理成 MAT。
%
% 固定的数据约定：
%   1. 离线算法输入 y(t) 优先来自原版 PX4 日志 topic：
%        vehicle_torque_setpoint.xyz  -> [Mx My Mz]
%        vehicle_thrust_setpoint.xyz  -> [Fx Fy Fz]
%      ControlAllocator::Run() 把这两个 3 维向量拼成 6 维 c[i]，然后调用：
%        _control_allocation[i]->setControlSetpoint(c[i])
%      所以这里保存的 control_sp 就是 allocator 真正收到的输入。
%
%   2. 如果日志里有多个 uORB instance：
%        vehicle_torque_setpoint instance 0 -> allocator matrix 0
%        vehicle_torque_setpoint instance 1 -> allocator matrix 1
%      main 里遇到实例数一致的 B 时会分别喂给每个 B。
%      如果日志实例数和待测 B 实例数不一致，只有在控制轴不重叠时才允许合并/拆分；
%      如果两个 instance 都能产生同一轴力/力矩，日志无法唯一说明该轴应该怎么重新分配。
%
%   3. 输出参考：
%        actuator_motors/servos 是原版 PX4 日志就有的最终执行器 topic。
%      parser 只把它们插值整理成 flightData.u_px4，用于离线结果对比。
%
%   4. 单位化不从日志输入判断。单位化只由待测 B 和算法配置决定：
%        control_allocation_scale = f(B, normalize_rpy)
%      输入 y(t) 保持日志里的 control_sp，不乘/除 scale。
%
% 输入变量：
%   LOG_PATH          .ulg 文件
%   COMMAND_INSTANCE  默认显示/窗口参考 instance，通常用 0
%   COMMAND_MAT       输出 MAT 文件
%   MODEL_NAME        机型名字，仅用于标记
%
% 输出：
%   flightData.control_sp          默认 instance 的 [Mx My Mz Fx Fy Fz]
%   flightData.setpoint_instances  原版输入 topic 的所有 instance
%   flightData.u_px4               PX4 在线执行器输出参考，按 motors 后 servos 拼接

fprintf('Parsing PX4 control-allocation replay data from:\n  %s\n', LOG_PATH);

%% 1. 临时导出目录
work_dir = tempname;
mkdir(work_dir);
cleanup_obj = onCleanup(@() rmdir(work_dir, 's')); %#ok<NASGU>

%% 2. 找 ulog2csv
ulog2csv_path = find_ulog2csv();
fprintf('Using ulog2csv:\n  %s\n', ulog2csv_path);

%% 3. 原版 PX4 输入：vehicle_torque_setpoint + vehicle_thrust_setpoint
torque_instances = read_xyz_topic_instances(ulog2csv_path, work_dir, LOG_PATH, ...
    'vehicle_torque_setpoint', false);
thrust_instances = read_xyz_topic_instances(ulog2csv_path, work_dir, LOG_PATH, ...
    'vehicle_thrust_setpoint', false);
setpoint_instances = build_control_setpoint_instances(torque_instances, thrust_instances);

%% 4. 选择默认 flightData，并固定 y/control_sp 语义
if ~isempty(setpoint_instances)
    selected_idx = find([setpoint_instances.instance] == COMMAND_INSTANCE, 1);

    if isempty(selected_idx)
        warning('setpoint instance %d not found. Use instance %d instead.', ...
            COMMAND_INSTANCE, setpoint_instances(1).instance);
        selected_idx = 1;
    end

    flightData = setpoint_instances(selected_idx);

else
    error(['No usable control-allocation input found. Need original PX4 topics ', ...
        'vehicle_torque_setpoint + vehicle_thrust_setpoint.']);
end

flightData.model = char(MODEL_NAME);
flightData.setpoint_instances = setpoint_instances;

%% 5. 可选输出参考：使用原版 actuator_motors/actuator_servos
[u_actuator, u_source] = read_actuator_outputs_if_available(ulog2csv_path, work_dir, LOG_PATH, ...
    flightData.timestamp_us, flightData.timestamp_sample_us);

flightData.u_px4 = u_actuator;
flightData.u_px4_source = u_source;
flightData.u_dim_px4 = size(flightData.u_px4, 2);

%% 6. 保存

if ~isfolder(fileparts(COMMAND_MAT))
    mkdir(fileparts(COMMAND_MAT));
end

save(COMMAND_MAT, 'flightData', '-v7.3');

fprintf('Saved command data:\n  %s\n', COMMAND_MAT);
fprintf('Input source: vehicle_torque_setpoint+vehicle_thrust_setpoint, selected instance: %d, samples: %d, dt_mean: %.6g s\n', ...
    flightData.instance, size(flightData.control_sp, 1), flightData.dt_mean_s);
fprintf('Setpoint instances: %s\n', mat2str(instances_to_ids(setpoint_instances)));
fprintf('PX4 output reference: %s, columns: %d\n', ...
    flightData.u_px4_source, size(flightData.u_px4, 2));

%% Local functions
function ulog2csv_path = find_ulog2csv()
% MATLAB 从 Finder/Dock 启动时 PATH 往往缺少 Python user bin。
% 所以这里按常见位置逐个找。
[status, out] = system('command -v ulog2csv');
ulog2csv_path = strtrim(out);

if status == 0 && ~isempty(ulog2csv_path)
    return;
end

home_dir = getenv('HOME');
candidate_paths = { ...
    fullfile(home_dir, 'Library', 'Python', '3.9', 'bin', 'ulog2csv'), ...
    '/Users/mch/Library/Python/3.9/bin/ulog2csv', ...
    '/opt/homebrew/bin/ulog2csv', ...
    '/usr/local/bin/ulog2csv'};

for i = 1:numel(candidate_paths)
    if isfile(candidate_paths{i})
        ulog2csv_path = candidate_paths{i};
        return;
    end
end

error('ulog2csv not found. Install pyulog, or edit candidate_paths in this parser.');
end

function csv_files = export_topic_csvs(ulog2csv_path, work_dir, log_path, topic, required)
% 导出一个 topic。required=false 时，原版日志缺自定义 topic 不报错。
cmd = sprintf('"%s" -m %s -o "%s" "%s"', ulog2csv_path, topic, work_dir, log_path);
[status, cmdout] = system(cmd);

if status ~= 0
    if required
        error('ulog2csv failed for %s:\n%s', topic, cmdout);
    end

    csv_files = [];
    return;
end

csv_files = dir(fullfile(work_dir, ['*_' topic '_*.csv']));

if isempty(csv_files)
    csv_files = dir(fullfile(work_dir, ['*_' topic '.csv']));
end

if isempty(csv_files) && required
    error('No %s CSV was exported from this log.', topic);
end
end

function instances = read_xyz_topic_instances(ulog2csv_path, work_dir, log_path, topic, required)
csv_files = export_topic_csvs(ulog2csv_path, work_dir, log_path, topic, required);

if isempty(csv_files)
    instances = repmat(empty_xyz_instance(topic), 1, 0);
    return;
end

instances = repmat(empty_xyz_instance(topic), 1, numel(csv_files));

for file_idx = 1:numel(csv_files)
    instance_id = instance_id_from_file(csv_files(file_idx).name, topic, file_idx);
    csv_path = fullfile(csv_files(file_idx).folder, csv_files(file_idx).name);
    instances(file_idx) = read_xyz_csv(csv_path, topic, instance_id);
end

[~, order] = sort([instances.instance]);
instances = instances(order);
end

function inst = read_xyz_csv(csv_path, topic, instance_id)
T = readtable(csv_path, 'VariableNamingRule', 'modify');
vars = T.Properties.VariableNames;
N = height(T);

inst = empty_xyz_instance(topic);
inst.instance = instance_id;
inst.timestamp_us = read_scalar_column(T, 'timestamp', nan(N, 1));
inst.timestamp_sample_us = read_scalar_column(T, 'timestamp_sample', inst.timestamp_us);
inst.dt_mean_s = mean_positive_dt(inst.timestamp_us);
inst.xyz = read_indexed_columns(T, vars, 'xyz', 3, nan(N, 3));
end

function setpoint_instances = build_control_setpoint_instances(torque_instances, thrust_instances)
% 复刻 ControlAllocator::Run() 里 c[i] 的构造：
%   c[i] = [torque.xyz thrust.xyz]
% allocator 只由 vehicle_torque_setpoint instance 0 的 callback 触发。
% 因此所有 allocator instance 都必须对齐到 instance 0 的运行时间轴：
%   matrix 0: 使用该时刻的 torque0/thrust0
%   matrix 1: copy 该时刻之前最新的 torque1/thrust1
if isempty(torque_instances)
    setpoint_instances = repmat(empty_setpoint_instance(), 1, 0);
    return;
end

reference_torque = find_matching_instance(torque_instances, 0);
ref_timestamp_us = reference_torque.timestamp_us;
ref_timestamp_sample_us = reference_torque.timestamp_sample_us;
N = numel(ref_timestamp_us);
setpoint_instances = repmat(empty_setpoint_instance(), 1, numel(torque_instances));

for i = 1:numel(torque_instances)
    torque = torque_instances(i);
    thrust = find_matching_instance(thrust_instances, torque.instance);

    if torque.instance == reference_torque.instance
        torque_xyz = fixed_width(torque.xyz, 3);
    else
        torque_xyz = align_rows_to_reference(torque.timestamp_us, torque.xyz, ref_timestamp_us, zeros(1, 3));
    end

    if isempty(thrust)
        thrust_xyz = zeros(N, 3);
    else
        thrust_xyz = align_rows_to_reference(thrust.timestamp_us, thrust.xyz, ref_timestamp_us, zeros(1, 3));
    end

    inst = empty_setpoint_instance();
    inst.instance = torque.instance;
    inst.timestamp_us = ref_timestamp_us;
    inst.timestamp_sample_us = ref_timestamp_sample_us;
    inst.dt_mean_s = reference_torque.dt_mean_s;
    inst.t = (0:(N - 1))' * inst.dt_mean_s;
    inst.control_sp = fixed_width_6([torque_xyz thrust_xyz]);
    setpoint_instances(i) = inst;
end
end

function [u_combined, source_name] = read_actuator_outputs_if_available(ulog2csv_path, work_dir, log_path, ...
    ref_timestamp_us, ref_timestamp_sample_us)
% 原版 PX4 输出 topic：
%   actuator_motors.control : motors，最多 12 列
%   actuator_servos.control : servos，最多 8 列
% 它们是最终发布到驱动层的执行器值，可能包含 slew/停桨/SITL actuator model。
motors = read_actuator_topic_if_available(ulog2csv_path, work_dir, log_path, ...
    'actuator_motors', 12, ref_timestamp_us, ref_timestamp_sample_us);
servos = read_actuator_topic_if_available(ulog2csv_path, work_dir, log_path, ...
    'actuator_servos', 8, ref_timestamp_us, ref_timestamp_sample_us);

u_combined = [motors.control servos.control];
u_combined = trim_trailing_all_nan_columns(u_combined);
source_name = 'actuator_motors+actuator_servos interpolated to command samples';
end

function out = read_actuator_topic_if_available(ulog2csv_path, work_dir, log_path, topic, control_count, ...
    ref_timestamp_us, ref_timestamp_sample_us)
out = struct('control', nan(numel(ref_timestamp_us), 0));
csv_files = export_topic_csvs(ulog2csv_path, work_dir, log_path, topic, false);

if isempty(csv_files)
    return;
end

csv_path = fullfile(csv_files(1).folder, csv_files(1).name);
T = readtable(csv_path, 'VariableNamingRule', 'modify');
vars = T.Properties.VariableNames;
N = height(T);
timestamp_us = read_scalar_column(T, 'timestamp', nan(N, 1));
timestamp_sample_us = read_scalar_column(T, 'timestamp_sample', timestamp_us);
control = read_indexed_columns(T, vars, 'control', control_count, nan(N, control_count));

source_time = choose_alignment_time(timestamp_us, timestamp_sample_us, ref_timestamp_sample_us);
ref_time = choose_reference_alignment_time(ref_timestamp_us, ref_timestamp_sample_us, timestamp_sample_us);
aligned = interpolate_rows_to_reference(source_time, control, ref_time);
aligned = trim_trailing_all_nan_columns(aligned);

out.control = aligned;
end

function values = read_indexed_columns(T, vars, prefix, count, default_values)
% 兼容 ulog2csv 常见列名：
%   y_0_ / y_0 / y0
values = default_values;

for idx = 0:(count - 1)
    candidates = {sprintf('%s_%d_', prefix, idx), sprintf('%s_%d', prefix, idx), sprintf('%s%d', prefix, idx)};
    col = find(ismember(vars, candidates), 1);

    if ~isempty(col)
        values(:, idx + 1) = T{:, col};
    end
end

values(~isfinite(values)) = default_values(~isfinite(values));
end

function values = read_scalar_column(T, name, default_values)
if ismember(name, T.Properties.VariableNames)
    values = T.(name);
else
    values = default_values;
end

values = double(values);
end

function aligned = align_rows_to_reference(source_time, source_values, ref_time, default_row)
% 输入 setpoint 对齐用 previous-hold，目的是复刻 Run() 中“保存最新 thrust”的行为。
% main benchmark 不再做 topic 对齐，它只消费这里已经整理好的 y(t)。
ref_time = double(ref_time(:));
source_time = double(source_time(:));
source_values = double(source_values);

if isempty(ref_time)
    aligned = zeros(0, size(source_values, 2));
    return;
end

if isempty(source_time) || isempty(source_values)
    aligned = repmat(default_row, numel(ref_time), 1);
    return;
end

aligned = repmat(default_row, numel(ref_time), 1);
source_idx = 1;

for i = 1:numel(ref_time)
    while source_idx < numel(source_time) && source_time(source_idx + 1) <= ref_time(i)
        source_idx = source_idx + 1;
    end

    if source_time(source_idx) <= ref_time(i) || source_idx == 1
        aligned(i, :) = source_values(source_idx, :);
    end
end
end

function aligned = interpolate_rows_to_reference(source_time, source_values, ref_time)
% 输出参考只用于和离线 u 画图/算差异；这里用线性插值对齐到 command 时间。
ref_time = double(ref_time(:));
source_time = double(source_time(:));
source_values = double(source_values);

if isempty(ref_time)
    aligned = zeros(0, size(source_values, 2));
    return;
end

if isempty(source_time) || isempty(source_values)
    aligned = nan(numel(ref_time), 0);
    return;
end

[source_time, unique_idx] = unique(source_time, 'stable');
source_values = source_values(unique_idx, :);
aligned = nan(numel(ref_time), size(source_values, 2));

for col = 1:size(source_values, 2)
    values = source_values(:, col);
    valid = isfinite(source_time) & isfinite(values);

    if nnz(valid) >= 2
        aligned(:, col) = interp1(source_time(valid), values(valid), ref_time, 'linear', nan);
    elseif nnz(valid) == 1
        aligned(:, col) = values(find(valid, 1));
    end
end
end

function t = choose_alignment_time(timestamp_us, timestamp_sample_us, ref_timestamp_sample_us)
if ~isempty(ref_timestamp_sample_us) && ~isempty(timestamp_sample_us) && any(isfinite(timestamp_sample_us))
    t = timestamp_sample_us;
else
    t = timestamp_us;
end
end

function t = choose_reference_alignment_time(ref_timestamp_us, ref_timestamp_sample_us, source_timestamp_sample_us)
if ~isempty(source_timestamp_sample_us) && any(isfinite(source_timestamp_sample_us)) ...
        && ~isempty(ref_timestamp_sample_us) && any(isfinite(ref_timestamp_sample_us))
    t = ref_timestamp_sample_us;
else
    t = ref_timestamp_us;
end
end

function dt_mean_s = mean_positive_dt(timestamp_us)
t_log = double(timestamp_us(:)) * 1e-6;

if numel(t_log) > 1
    dt_values = diff(t_log);
    dt_values = dt_values(isfinite(dt_values) & dt_values > 0);

    if isempty(dt_values)
        dt_mean_s = 0.01;
    else
        dt_mean_s = mean(dt_values);
    end
else
    dt_mean_s = 0.01;
end
end

function inst = find_matching_instance(instances, instance_id)
inst = [];

if isempty(instances)
    return;
end

idx = find([instances.instance] == instance_id, 1);

if isempty(idx)
    idx = find([instances.instance] == 0, 1);
end

if isempty(idx)
    idx = 1;
end

inst = instances(idx);
end

function instance_id = instance_id_from_file(file_name, topic, file_idx)
token = regexp(file_name, ['_' topic '_(\d+)\.csv$'], 'tokens', 'once');

if isempty(token)
    instance_id = file_idx - 1;
else
    instance_id = str2double(token{1});
end
end

function ids = instances_to_ids(instances)
if isempty(instances)
    ids = [];
else
    ids = [instances.instance];
end
end

function values = fixed_width(values, width)
values = double(values);
values(~isfinite(values)) = 0;

if size(values, 2) < width
    values(:, end+1:width) = 0;
elseif size(values, 2) > width
    values = values(:, 1:width);
end
end

function Y = fixed_width_6(Y)
Y = fixed_width(Y, 6);
end

function values = trim_trailing_all_nan_columns(values)
if isempty(values)
    return;
end

last_col = find(any(isfinite(values), 1), 1, 'last');

if isempty(last_col)
    values = zeros(size(values, 1), 0);
else
    values = values(:, 1:last_col);
end
end

function inst = empty_xyz_instance(~)
inst = struct();
inst.instance = 0;
inst.timestamp_us = [];
inst.timestamp_sample_us = [];
inst.dt_mean_s = nan;
inst.xyz = [];
end

function inst = empty_setpoint_instance()
inst = struct();
inst.instance = 0;
inst.timestamp_us = [];
inst.timestamp_sample_us = [];
inst.t = [];
inst.dt_mean_s = nan;
inst.control_sp = [];
end
