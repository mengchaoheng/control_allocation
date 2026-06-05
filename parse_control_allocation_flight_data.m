%% Parse ULog topics needed by main_control_allocation_benchmark.m
%
% 输入变量来自 main：
%   LOG_PATH       : .ulg 文件完整路径
%   COMMAND_MAT    : 输出 MAT 路径
%   PYULOG_BIN_DIR : 可选，ulog2csv/ulog_info/ulog_params 所在目录；空则从 PATH 自动找
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
%   airframe_id            : SYS_AUTOSTART，例如 22005
%   log_model_name         : string(airframe_id)
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

if ~exist('PYULOG_BIN_DIR', 'var')
    PYULOG_BIN_DIR = "";
end

ulog2csv_path = find_ulog_tool('ulog2csv', PYULOG_BIN_DIR);
topics = { ...
    'vehicle_torque_setpoint', ...
    'vehicle_thrust_setpoint', ...
    'actuator_motors', ...
    'actuator_servos', ...
    'allocation_value'};

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
[airframe_id, log_model_name] = read_airframe_id(LOG_PATH, PYULOG_BIN_DIR);
fprintf('Airframe: %s\n', log_model_name);

save_names = [save_names, {'v_sp_t', 'log_v_sp_instances', 'u_px4', 'u_px4_source', ...
    'airframe_id', 'log_model_name'}];

save(COMMAND_MAT, save_names{:}, '-v7.3');

fprintf('Saved raw topic MAT:\n  %s\n', COMMAND_MAT);
fprintf('Variables: %s\n', strjoin(save_names, ', '));

function [airframe_id, log_model_name] = read_airframe_id(log_path, pyulog_bin_dir)
airframe_id = nan;

ulog_info_path = find_ulog_tool('ulog_info', pyulog_bin_dir);
[status, out] = system(sprintf('"%s" "%s"', ulog_info_path, log_path));

if status == 0
    token = regexp(out, 'Airframe:\s*([0-9]+)', 'tokens', 'once');

    if ~isempty(token)
        airframe_id = str2double(token{1});
    end
end

if ~isfinite(airframe_id)
    ulog_params_path = find_ulog_tool('ulog_params', pyulog_bin_dir);
    [status, out] = system(sprintf('"%s" -i -f csv "%s"', ulog_params_path, log_path));

    if status == 0
        token = regexp(out, '(?m)^SYS_AUTOSTART,([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once');

        if ~isempty(token)
            airframe_id = str2double(token{1});
        end
    end
end

if isfinite(airframe_id)
    log_model_name = string(round(airframe_id));
else
    log_model_name = "unknown";
end
end

function tool_path = find_ulog_tool(tool_name, pyulog_bin_dir)
pyulog_bin_dir = string(pyulog_bin_dir);

if strlength(pyulog_bin_dir) > 0
    tool_path = fullfile(char(pyulog_bin_dir), tool_name);

    if isfile(tool_path)
        return;
    end

    error('%s not found in PYULOG_BIN_DIR: %s', tool_name, pyulog_bin_dir);
end

[status, out] = system(['command -v ' tool_name]);
tool_path = strtrim(out);

if status ~= 0 || isempty(tool_path)
    error('%s not found. Set PYULOG_BIN_DIR in main_control_allocation_benchmark.m or add it to PATH.', tool_name);
end

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
