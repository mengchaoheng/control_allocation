%% Parse ULog topics needed by main_control_allocation_benchmark.m
%
% 输入变量来自 main：
%   LOG_PATH    : .ulg 文件完整路径
%   COMMAND_MAT : 输出 MAT 路径
%
% 输出 MAT 保存原始 PX4 topic 表，变量名和日志 topic 名保持一致：
%   vehicle_torque_setpoint_0
%   vehicle_thrust_setpoint_0
%   vehicle_torque_setpoint_1      若日志存在
%   vehicle_thrust_setpoint_1      若日志存在
%   actuator_motors_0              若日志存在
%   actuator_servos_0              若日志存在
%
% 同时直接构造 main 使用的变量：
%   v_sp_t                 : vehicle_torque_setpoint_0.timestamp, seconds
%   log_v_sp_instances{i}  : 6 x N, [Mx My Mz Fx Fy Fz]
%   u_px4                  : actuator_motors_0 + actuator_servos_0, 插值到 v_sp_t
%
% v_sp 构造关系：
%   v(1:3,:) = vehicle_torque_setpoint_*.xyz
%   v(4:6,:) = vehicle_thrust_setpoint_*.xyz

fprintf('Parsing ULog topics for control allocation:\n  %s\n', LOG_PATH);

if ~isfolder(fileparts(COMMAND_MAT))
    mkdir(fileparts(COMMAND_MAT));
end

work_dir = tempname;
mkdir(work_dir);
cleanup_obj = onCleanup(@() rmdir(work_dir, 's')); %#ok<NASGU>

ulog2csv_path = find_ulog2csv();
topics = { ...
    'vehicle_torque_setpoint', ...
    'vehicle_thrust_setpoint', ...
    'actuator_motors', ...
    'actuator_servos'};

save_names = {};

for topic_idx = 1:numel(topics)
    topic = topics{topic_idx};
    cmd = sprintf('"%s" -m %s -o "%s" "%s"', ulog2csv_path, topic, work_dir, LOG_PATH);
    fprintf('  %s\n', cmd);
    [status, cmdout] = system(cmd);

    if status ~= 0
        fprintf('  skip %s: %s\n', topic, strtrim(cmdout));
        continue;
    end

    csv_files = dir(fullfile(work_dir, ['*_' topic '_*.csv']));
    if isempty(csv_files)
        csv_files = dir(fullfile(work_dir, ['*_' topic '.csv']));
    end

    for file_idx = 1:numel(csv_files)
        csv_path = fullfile(csv_files(file_idx).folder, csv_files(file_idx).name);
        instance_id = instance_id_from_csv_name(csv_files(file_idx).name, topic);
        var_name = sprintf('%s_%d', topic, instance_id);
        T = readtable(csv_path, 'VariableNamingRule', 'modify'); %#ok<NASGU>
        eval([var_name ' = T;']);
        save_names{end+1} = var_name; %#ok<SAGROW>
        fprintf('    saved variable: %s (%d rows)\n', var_name, height(T));
    end
end

save_names = unique(save_names, 'stable');

v_sp_t = double(vehicle_torque_setpoint_0.timestamp(:)) * 1e-6;

vehicle_torque_setpoint_0_xyz = [
    vehicle_torque_setpoint_0.xyz_0_(:)';
    vehicle_torque_setpoint_0.xyz_1_(:)';
    vehicle_torque_setpoint_0.xyz_2_(:)'];

vehicle_thrust_setpoint_0_t = double(vehicle_thrust_setpoint_0.timestamp(:)) * 1e-6;
vehicle_thrust_setpoint_0_xyz_raw = [
    vehicle_thrust_setpoint_0.xyz_0_(:), ...
    vehicle_thrust_setpoint_0.xyz_1_(:), ...
    vehicle_thrust_setpoint_0.xyz_2_(:)];
vehicle_thrust_setpoint_0_xyz = interp1(vehicle_thrust_setpoint_0_t, ...
    vehicle_thrust_setpoint_0_xyz_raw, v_sp_t, 'previous', 'extrap')';

log_v_sp_instances = {[vehicle_torque_setpoint_0_xyz; vehicle_thrust_setpoint_0_xyz]};

if exist('vehicle_torque_setpoint_1', 'var')
    vehicle_torque_setpoint_1_t = double(vehicle_torque_setpoint_1.timestamp(:)) * 1e-6;
    vehicle_torque_setpoint_1_xyz_raw = [
        vehicle_torque_setpoint_1.xyz_0_(:), ...
        vehicle_torque_setpoint_1.xyz_1_(:), ...
        vehicle_torque_setpoint_1.xyz_2_(:)];
    vehicle_torque_setpoint_1_xyz = interp1(vehicle_torque_setpoint_1_t, ...
        vehicle_torque_setpoint_1_xyz_raw, v_sp_t, 'previous', 'extrap')';

    if exist('vehicle_thrust_setpoint_1', 'var')
        vehicle_thrust_setpoint_1_t = double(vehicle_thrust_setpoint_1.timestamp(:)) * 1e-6;
        vehicle_thrust_setpoint_1_xyz_raw = [
            vehicle_thrust_setpoint_1.xyz_0_(:), ...
            vehicle_thrust_setpoint_1.xyz_1_(:), ...
            vehicle_thrust_setpoint_1.xyz_2_(:)];
        vehicle_thrust_setpoint_1_xyz = interp1(vehicle_thrust_setpoint_1_t, ...
            vehicle_thrust_setpoint_1_xyz_raw, v_sp_t, 'previous', 'extrap')';
    else
        vehicle_thrust_setpoint_1_xyz = zeros(3, numel(v_sp_t));
    end

    log_v_sp_instances{end+1} = [vehicle_torque_setpoint_1_xyz; vehicle_thrust_setpoint_1_xyz];
end

u_px4 = [];
if exist('actuator_motors_0', 'var')
    actuator_motors_0_t = double(actuator_motors_0.timestamp(:)) * 1e-6;
    actuator_motors_0_cols = startsWith(actuator_motors_0.Properties.VariableNames, 'control_');
    actuator_motors_0_control = actuator_motors_0{:, actuator_motors_0_cols};
    u_px4 = [u_px4, interp1(actuator_motors_0_t, actuator_motors_0_control, v_sp_t, 'previous', 'extrap')];
end

if exist('actuator_servos_0', 'var')
    actuator_servos_0_t = double(actuator_servos_0.timestamp(:)) * 1e-6;
    actuator_servos_0_cols = startsWith(actuator_servos_0.Properties.VariableNames, 'control_');
    actuator_servos_0_control = actuator_servos_0{:, actuator_servos_0_cols};
    u_px4 = [u_px4, interp1(actuator_servos_0_t, actuator_servos_0_control, v_sp_t, 'previous', 'extrap')];
end

if isempty(u_px4)
    u_px4 = nan(numel(v_sp_t), 0);
else
    u_px4 = u_px4(:, any(isfinite(u_px4), 1));
end

u_px4_source = 'actuator_motors_0 + actuator_servos_0';
save_names = [save_names, {'v_sp_t', 'log_v_sp_instances', 'u_px4', 'u_px4_source'}];

save(COMMAND_MAT, save_names{:}, '-v7.3');

fprintf('Saved raw topic MAT:\n  %s\n', COMMAND_MAT);
fprintf('Variables: %s\n', strjoin(save_names, ', '));

function ulog2csv_path = find_ulog2csv()
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

error('ulog2csv not found.');
end

function instance_id = instance_id_from_csv_name(file_name, topic)
pattern = ['_' topic '_(\d+)\.csv$'];
token = regexp(file_name, pattern, 'tokens', 'once');

if isempty(token)
    instance_id = 0;
else
    instance_id = str2double(token{1});
end
end
