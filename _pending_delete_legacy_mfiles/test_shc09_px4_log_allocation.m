clearvars;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpath(script_dir));

mat_path = '/Users/mch/Proj/control_allocation/15_58_23.mat';
aligned_cache_dir = fullfile(script_dir, 'cache');
[~, mat_name] = fileparts(mat_path);
aligned_data_path = fullfile(aligned_cache_dir, [mat_name '_aligned_vu.mat']);
legacy_aligned_data_path = fullfile(aligned_cache_dir, [mat_name '_aligned_yu.mat']);
plot_cache_path = fullfile(aligned_cache_dir, 'last_shc09_allocation_plot_data.mat');
rebuild_aligned_data = false; % true 时才重新解析原始 ulog MAT 并覆盖 aligned_data_path。
time_window_s = [];      % [] = 全部数据；否则只跑这个时间窗口，例如 [10 20]。
use_restoring = true;
export_cpp_replay_inputs_enabled = true;
compare_cpp_replay_outputs = true;
auto_run_cpp_replay = false;
cpp_replay_binary = fullfile(script_dir, 'alloc_cpp', 'build', 'shc09_log_replay');
cpp_replay_input_dir = fullfile(script_dir, 'results', 'cpp_inputs');
cpp_replay_output_dir = fullfile(script_dir, 'results', 'cpp_outputs');
cpp_results = struct();

allocation_method_selection = { ...
    'px4_inv', ...
    'inv_Bnorm', ...
    'inv_Bpar', ...
    'split_inv', ...
    'pca_dir_bpar', ...
    'pca_dpscaled_bpar', ...
    'split_pca_dir', ...
    'split_pca_dpscaled'};
% allocation_method_selection 可设为 'all'、数字索引，或方法名 cell。
% 新增方法时，在 get_allocation_method_catalog 和下面的 results 生成段各加一处。
reports_to_run = {'allocator_summary', 'method_diff'};
plot_reference_name = 'inv_Bpar';
tie_opts = struct(); % 不再做 original/tiebreak 对照，DP/PCA 走默认 simplex 后端。

allocation_method_catalog = get_allocation_method_catalog();
allocation_methods_to_run = resolve_allocation_method_selection( ...
    allocation_method_selection, allocation_method_catalog);
method_enabled = @(name) any(strcmp(allocation_methods_to_run, name));

% 本脚本主线：固定 SHC09 的 B，从日志取输入 v，再切换不同分配策略查看输出 u。
% SHC09: v = B*u, v=[Mx My Mz Fx Fy Fz]', u=[motor0 servo0 ... servo5]'.
% 对分配输出 u 来说：
%   同一行内各执行器的相对值更直接决定“这个控制轴由谁分担”；
%   同一列内力/力矩的相对值决定“某个执行器带来的跨轴耦合和补偿”。
axis_labels = {'Mx', 'My', 'Mz', 'Fx', 'Fy', 'Fz'};
channel_labels = {'motor0', 'servo0', 'servo1', 'servo2', 'servo3', 'servo4', 'servo5'};
B = [
    0        -0.5010  -0.2505   0.2505   0.5010   0.2505  -0.2505
    0         0        0.4339   0.4339   0       -0.4339  -0.4339
    0.3250    0.2070   0.2070   0.2070   0.2070   0.2070   0.2070
    0         0        0        0        0        0        0
    0         0        0        0        0        0        0
   -6.5000    0        0        0        0        0        0
];
umin = [0; -1; -1; -1; -1; -1; -1];
umax = [1;  1;  1;  1;  1;  1;  1];

[v, u_reference, t, cache_status] = load_or_make_aligned_vu( ...
    mat_path, aligned_data_path, legacy_aligned_data_path, rebuild_aligned_data);
u_reference_label = 'u_log';

t = t(:)';
if isempty(time_window_s)
    time_ok = true(size(t));
else
    if numel(time_window_s) ~= 2
        error('test_shc09:BadTimeWindow', 'time_window_s must be [] or [t_start t_end].');
    end
    time_ok = t >= time_window_s(1) & t <= time_window_s(2);
end
v = v(:, time_ok);
u_reference = u_reference(:, time_ok);
t = t(time_ok);
if isempty(t)
    error('test_shc09:EmptyTimeWindow', 'No samples found in time_window_s.');
end
sample_count = size(v, 2);
zero_errout = zeros(1, sample_count);
unit_scale = ones(1, sample_count);

% PX4 mix normalization: mix_norm = B^+ * D^-1, 等价用 B_norm = D*B做分配。
mix_px4_raw = px4_geninv(B);
scale = px4_shc09_scale_from_mix(mix_px4_raw);
mix_px4_norm = mix_px4_raw ./ scale';
B_norm = diag(scale) * B;

active_rows = find(any(abs(B_norm) > 1e-10, 2));
B_par = B_norm(active_rows, :);
v_par = v(active_rows, :);
y_par = v_par; % 兼容旧 plot cache 字段名。
axis_par = axis_labels(active_rows);

results = struct();

% 1. PX4 日志对应方法：原始 B 取 geninv 后按 PX4 scale 归一化 mix。
% 这里用于验证日志中的 PX4 inv 输出，不作为后续 PCA 对比的 B。
if method_enabled('px4_inv')
    u_px4_raw = mix_px4_norm * v;
    u_px4 = clamp_cols(u_px4_raw, umin, umax);
    results.px4_inv = make_result(u_px4, zero_errout, unit_scale);
end

% 2. 不拆力/力矩：完整 B_norm  
if method_enabled('inv_Bnorm')
    mix_inv_Bnorm = pinv(B_norm);
    u_inv_Bnorm_raw = mix_inv_Bnorm * v;
    u_inv_Bnorm = clamp_cols(u_inv_Bnorm_raw, umin, umax);
    results.inv_Bnorm = make_result(u_inv_Bnorm, zero_errout, unit_scale);
end

% 3. 删除零行后的 B_par。B_par 只是删掉 Fx/Fy 两个全零行。
% 3.1 伪逆法。inv_Bnorm 和 inv_Bpar 应当几乎一致；
mix_inv_Bpar = pinv(B_par);
u_inv_Bpar_raw = mix_inv_Bpar * v_par;
u_inv_Bpar = clamp_cols(u_inv_Bpar_raw, umin, umax);
if method_enabled('inv_Bpar')
    results.inv_Bpar = make_result(u_inv_Bpar, zero_errout, unit_scale);
end

if isempty(u_reference)
    u_reference = u_inv_Bpar;
end

% 4. split：先分配力，再补偿力执行器附带的力矩，最后分配剩余力矩。
%
% 这个分解依赖 SHC09 的结构：
%   v = B*u, v=[tau; f], u=[u_force; u_servo]
%   tau = B_torque_force*u_force + B_torque_servo*u_servo
%   f   = B_force_force*u_force + 0*u_servo
%
% 也就是舵面 servo 不产生力，B_force_servo=0。因此可以先用力通道算
% u_force。力通道可以是一对一除法，也可以是 pinv(B_force)*f 后 clamp；
% 如果发生限幅，后续补偿必须使用 clamp 后的真实 u_force，而不是 raw。
%
% 接着计算这个真实 u_force 已经带来的力矩：
%   tau_from_force = B_torque_force*u_force
% 并从期望力矩中扣除：
%   tau_left = tau_cmd - tau_from_force
% 最后只用不产生力的 servo 分配 tau_left。
%
% 这样做在本构型下相当于固定维度的结构化解耦。它避免把 Fz 和力矩
% 一起丢进 DP/PCA 后，由 motor0 的 Fz-Mz 耦合引入不必要的 LP 退化。
% 若未来舵面也进入力模型，即 B_force_servo ~= 0，则该分解不再严格成立。
fz_row = find(strcmp(axis_par, 'Fz'));
torque_rows = find(ismember(axis_par, {'Mx', 'My', 'Mz'}));
motor_col = find(abs(B_par(fz_row, :)) > 1e-10);
servo_cols = setdiff(1:numel(umin), motor_col, 'stable');

u_force = zeros(numel(umin), sample_count);
u_force(motor_col, :) = v_par(fz_row, :) / B_par(fz_row, motor_col);
u_force = clamp_cols(u_force, umin, umax);

% motor0 已经用于 Fz；剩余力矩只交给 6 个 servo。
B_torque_servo = B_par(torque_rows, servo_cols);
umin_servo = umin(servo_cols);
umax_servo = umax(servo_cols);

% Fz 的 motor0 输出也会带来 Mz，所以先从三维力矩期望里扣除这部分。
v_torque_from_motor = B_par(torque_rows, :) * u_force;
v_torque_left = v_par(torque_rows, :) - v_torque_from_motor;

% 4.1 split / inv：只对剩余三维力矩做 3x6 pinv 分配。
if method_enabled('split_inv')
    mix_torque_inv = pinv(B_torque_servo);
    u_torque_inv_raw = mix_torque_inv * v_torque_left;
    u_torque_inv = clamp_cols(u_torque_inv_raw, umin_servo, umax_servo);
    u_split_inv = add_servo_solution(u_force, servo_cols, u_torque_inv);
    results.split_inv = make_result(u_split_inv, zero_errout, unit_scale);
end

% 5. PCA/DP 方法。不再把 simplex 后端作为测试维度；方法名体现 B_par / split。
if method_enabled('pca_dir_bpar')
    results.pca_dir_bpar = run_dp_series( ...
        v_par, B_par, umin, umax, use_restoring, tie_opts);
end

if method_enabled('pca_dpscaled_bpar')
    results.pca_dpscaled_bpar = run_dpscaled_series( ...
        v_par, B_par, umin, umax, use_restoring, tie_opts);
end

if method_enabled('split_pca_dir')
    torque_dir = run_dp_series( ...
        v_torque_left, B_torque_servo, umin_servo, umax_servo, ...
        use_restoring, tie_opts);
    results.split_pca_dir = make_result( ...
        add_servo_solution(u_force, servo_cols, torque_dir.u), ...
        torque_dir.errout, torque_dir.scale);
end

if method_enabled('split_pca_dpscaled')
    torque_scaled = run_dpscaled_series( ...
        v_torque_left, B_torque_servo, umin_servo, umax_servo, ...
        use_restoring, tie_opts);
    results.split_pca_dpscaled = make_result( ...
        add_servo_solution(u_force, servo_cols, torque_scaled.u), ...
        torque_scaled.errout, torque_scaled.scale);
end

fprintf('\nB_norm rows kept for B_par: [%s] = [%s]\n', join_num(active_rows), strjoin(axis_par, ' '));
fprintf('aligned v/u cache: %s\n', cache_status);
fprintf('u reference for rms/black line: %s\n', u_reference_label);
if isempty(time_window_s)
    fprintf('time window: all data %.6g to %.6g s, samples=%d\n', ...
        t(1), t(end), sample_count);
else
    fprintf('time window: %.6g to %.6g s, samples=%d\n', ...
        time_window_s(1), time_window_s(2), sample_count);
end
fprintf('allocation methods: %s\n', strjoin(allocation_methods_to_run, ', '));
fprintf('restoring: %s\n', string(use_restoring));

if isempty(fieldnames(results))
    error('test_shc09:NoMethodsRun', 'No allocation methods were run.');
end

if ~isfield(results, plot_reference_name)
    plot_reference_name = first_result_name(results);
    fprintf('plot/diff reference changed to first enabled method: %s\n', plot_reference_name);
end

if has_report(reports_to_run, 'allocator_summary')
    print_summary(results, B_par, v_par, u_reference, umin, umax, u_reference_label);
end

if has_report(reports_to_run, 'method_diff')
    print_difference_summary(results, plot_reference_name, B_par, v_par, umin, umax);
end

if export_cpp_replay_inputs_enabled
    export_cpp_replay_inputs(cpp_replay_input_dir, ...
        B_par, v_par, umin, umax, ...
        B_torque_servo, v_torque_left, u_force, umin_servo, umax_servo, servo_cols, t);
    fprintf('C++ replay input dir: %s\n', cpp_replay_input_dir);
end

if compare_cpp_replay_outputs
    if auto_run_cpp_replay
        if ~run_cpp_replay_binary(cpp_replay_binary, script_dir)
            fprintf('C++ replay comparison skipped because the binary did not run.\n');
            cpp_results = struct();
        else
            cpp_results = load_cpp_replay_results(cpp_replay_output_dir);
        end
    else
        cpp_results = load_cpp_replay_results(cpp_replay_output_dir);
    end

    if isempty(fieldnames(cpp_results))
        fprintf('C++ replay outputs not found in: %s\n', cpp_replay_output_dir);
        fprintf('From alloc_cpp/build run: ./shc09_log_replay\n');
    else
        loaded_cpp_methods = fieldnames(cpp_results);
        fprintf('C++ replay output dir: %s\n', cpp_replay_output_dir);
        fprintf('C++ replay loaded methods: %s\n', strjoin(loaded_cpp_methods.', ', '));
        report_cpp_replay_differences(results, cpp_results, B_par, v_par);
    end
end

save(plot_cache_path, ...
    't', 'u_reference', 'u_reference_label', 'results', ...
    'channel_labels', 'B_par', 'v_par', 'y_par', 'axis_par', 'plot_reference_name', ...
    'allocation_method_selection', 'allocation_method_catalog', ...
    'allocation_methods_to_run', 'reports_to_run', 'use_restoring', 'tie_opts', ...
    'aligned_data_path', 'export_cpp_replay_inputs_enabled', ...
    'compare_cpp_replay_outputs', 'auto_run_cpp_replay', 'cpp_replay_binary', ...
    'cpp_replay_input_dir', 'cpp_replay_output_dir', 'cpp_results');
fprintf('plot cache: %s\n', plot_cache_path);

% test 跑完立即画图；之后也可以直接运行 plot_all_results 重画上次结果。
plot_all_results(t, u_reference, results, channel_labels, B_par, v_par, axis_par, plot_reference_name, u_reference_label, cpp_results);

function out = make_result(u, errout, scale)
    out = struct('u', u, 'errout', errout, 'scale', scale);
end

function catalog = get_allocation_method_catalog()
    catalog = { ...
        'px4_inv', ...            % PX4 geninv + PX4 scale，验证日志 mixer 风格输出。
        'inv_Bnorm', ...          % pinv(B_norm)*v + clamp。
        'inv_Bpar', ...           % pinv(B_par)*v_par + clamp，作为 B_par 常用参考。
        'split_inv', ...          % Fz 直接分配给 motor0，剩余 torque 用 pinv 分配给 servo。
        'pca_dir_bpar', ...       % B_par 上的 DP_LPCA。
        'pca_dpscaled_bpar', ...  % B_par 上的 DPscaled_LPCA。
        'split_pca_dir', ...      % split 后剩余 torque 用 DP_LPCA。
        'split_pca_dpscaled'};    % split 后剩余 torque 用 DPscaled_LPCA。
end

function methods = resolve_allocation_method_selection(selection, catalog)
    if isempty(selection)
        methods = catalog;
        return;
    end

    if ischar(selection) || isstring(selection)
        selection = cellstr(selection);

        if numel(selection) == 1 && strcmpi(selection{1}, 'all')
            methods = catalog;
            return;
        end

        methods = validate_allocation_methods(selection, catalog);
        return;
    end

    if isnumeric(selection)
        methods = catalog(selection);
        return;
    end

    if iscell(selection)
        methods = validate_allocation_methods(selection, catalog);
        return;
    end

    error('test_shc09:BadAllocatorSelection', ...
          'allocation_method_selection must be ''all'', numeric indices, or method names.');
end

function methods = validate_allocation_methods(methods, catalog)
    canonical = cell(1, numel(methods));

    for i = 1:numel(methods)
        idx = find(strcmpi(catalog, char(methods{i})), 1);
        if isempty(idx)
            error('test_shc09:UnknownAllocatorMethod', ...
                  'Unknown allocator method "%s". Use one of: %s', ...
                  char(methods{i}), strjoin(catalog, ', '));
        end
        canonical{i} = catalog{idx};
    end

    [~, unique_idx] = unique(canonical, 'stable');
    methods = canonical(sort(unique_idx));
end

function enabled = has_report(reports_to_run, report_name)
    if ischar(reports_to_run) || isstring(reports_to_run)
        reports_to_run = cellstr(reports_to_run);
    end

    enabled = any(strcmpi(reports_to_run, 'all')) || any(strcmpi(reports_to_run, report_name));
end

function name = first_result_name(results)
    names = fieldnames(results);
    name = names{1};
end

function export_cpp_replay_inputs(input_dir, B_par, v_par, umin, umax, B_torque_servo, v_torque_left, u_force, umin_servo, umax_servo, servo_cols, t)
    if exist(input_dir, 'dir') ~= 7
        mkdir(input_dir);
    end

    writematrix(B_par, fullfile(input_dir, 'shc09_log_bpar.csv'));
    writematrix(v_par', fullfile(input_dir, 'shc09_log_v_par.csv'));
    writematrix(umin(:)', fullfile(input_dir, 'shc09_log_umin.csv'));
    writematrix(umax(:)', fullfile(input_dir, 'shc09_log_umax.csv'));
    writematrix(B_torque_servo, fullfile(input_dir, 'shc09_log_b_torque_servo.csv'));
    writematrix(v_torque_left', fullfile(input_dir, 'shc09_log_v_torque_left.csv'));
    writematrix(u_force', fullfile(input_dir, 'shc09_log_u_force.csv'));
    writematrix(umin_servo(:)', fullfile(input_dir, 'shc09_log_umin_servo.csv'));
    writematrix(umax_servo(:)', fullfile(input_dir, 'shc09_log_umax_servo.csv'));
    writematrix(servo_cols(:)' - 1, fullfile(input_dir, 'shc09_log_servo_cols_zero_based.csv'));
    writematrix(t(:), fullfile(input_dir, 'shc09_log_t.csv'));
end

function ok = run_cpp_replay_binary(cpp_binary, script_dir)
    ok = false;
    if exist(cpp_binary, 'file') ~= 2
        warning('test_shc09:MissingCppReplay', ...
            'C++ replay binary not found: %s. Build from alloc_cpp/build with: cmake .. && make shc09_log_replay', cpp_binary);
        return;
    end

    cmd = sprintf('"%s" "%s"', cpp_binary, script_dir);
    [status, output] = system(cmd);
    if status ~= 0
        warning('test_shc09:CppReplayFailed', 'C++ replay failed:\n%s', output);
        return;
    end

    output = strtrim(output);
    if ~isempty(output)
        fprintf('\n%s\n', output);
    end
    ok = true;
end

function cpp_results = load_cpp_replay_results(output_dir)
    methods = {'pca_dir_bpar', 'pca_dpscaled_bpar', 'split_pca_dir', 'split_pca_dpscaled'};
    cpp_results = struct();

    for i = 1:numel(methods)
        method = methods{i};
        u_file = fullfile(output_dir, ['output_cpp_shc09_log_' method '.csv']);
        err_file = fullfile(output_dir, ['output_cpp_shc09_log_' method '_errout.csv']);
        rho_file = fullfile(output_dir, ['output_cpp_shc09_log_' method '_rho.csv']);

        if exist(u_file, 'file') ~= 2
            continue;
        end

        u = readmatrix(u_file)';
        errout = zeros(1, size(u, 2));
        scale = ones(1, size(u, 2));
        if exist(err_file, 'file') == 2
            errout = readmatrix(err_file)';
            errout = errout(:)';
        end
        if exist(rho_file, 'file') == 2
            scale = readmatrix(rho_file)';
            scale = scale(:)';
        end
        cpp_results.(method) = make_result(u, errout, scale);
    end
end

function report_cpp_replay_differences(matlab_results, cpp_results, B_eval, v_eval)
    methods = {'pca_dir_bpar', 'pca_dpscaled_bpar', 'split_pca_dir', 'split_pca_dpscaled'};

    fprintf('\nC++ replay differences vs MATLAB\n');
    fprintf('%-22s %10s %12s %12s %12s %12s %10s %10s\n', ...
        'method', 'ok_common', 'rms_du_ok', 'max_du_ok', 'max_dBu_ok', 'max_du_all', 'cpp_err', 'mat_err');

    for i = 1:numel(methods)
        method = methods{i};
        if ~isfield(matlab_results, method)
            fprintf('%-22s skipped: MATLAB method not enabled\n', method);
            continue;
        end
        if ~isfield(cpp_results, method)
            fprintf('%-22s skipped: C++ output missing\n', method);
            continue;
        end

        u_mat = matlab_results.(method).u;
        u_cpp = cpp_results.(method).u;
        n = min([size(u_mat, 2), size(u_cpp, 2), size(v_eval, 2)]);
        good = all(isfinite(u_mat(:, 1:n)), 1) & all(isfinite(u_cpp(:, 1:n)), 1);

        if ~any(good)
            fprintf('%-22s no finite common samples\n', method);
            continue;
        end

        du = u_cpp(:, good) - u_mat(:, good);
        mat_errout = matlab_results.(method).errout(1:n);
        cpp_errout = cpp_results.(method).errout(1:n);
        common_ok = good & (mat_errout == 0) & (cpp_errout == 0);

        if any(common_ok)
            du_ok = u_cpp(:, common_ok) - u_mat(:, common_ok);
            dBu_ok = B_eval * u_cpp(:, common_ok) - B_eval * u_mat(:, common_ok);
            rms_du_ok = rms_vec(du_ok);
            max_du_ok = max(abs(du_ok), [], 'all');
            max_dBu_ok = max(abs(dBu_ok), [], 'all');
        else
            rms_du_ok = NaN;
            max_du_ok = NaN;
            max_dBu_ok = NaN;
        end

        fprintf('%-22s %10d %12.4g %12.4g %12.4g %12.4g %10d %10d\n', ...
            method, nnz(common_ok), rms_du_ok, max_du_ok, max_dBu_ok, ...
            max(abs(du), [], 'all'), nnz(cpp_errout ~= 0), nnz(mat_errout ~= 0));
    end
end

function [v, u_log, t, cache_status] = load_or_make_aligned_vu(mat_path, aligned_data_path, legacy_aligned_data_path, rebuild_aligned_data)
    % 控制分配测试只依赖对齐后的 v/u_log/t。
    % 默认只读 aligned_data_path；设 rebuild_aligned_data=true 才重新解析原始 ulog MAT。
    extraction_version = 2;

    if nargin < 4
        rebuild_aligned_data = false;
    end

    if ~rebuild_aligned_data
        load_path = find_existing_aligned_data(aligned_data_path, legacy_aligned_data_path);
        if ~isempty(load_path)
            [v, u_log, t, cache_info] = load_aligned_vu_cache(load_path);
            cache_status = sprintf('loaded %s', load_path);

            if ~strcmp(load_path, aligned_data_path)
                save_aligned_vu_cache(aligned_data_path, v, u_log, t, cache_info);
                cache_status = sprintf('%s; migrated %s', cache_status, aligned_data_path);
            end
            return;
        end
    end

    source_file = dir(mat_path);
    if isempty(source_file)
        error('找不到对齐数据缓存，也找不到原始日志 MAT 文件: %s', mat_path);
    end

    [v, u_log, t] = load_log_vu(mat_path);

    cache_info = struct();
    cache_info.source_path = mat_path;
    cache_info.source_bytes = source_file.bytes;
    cache_info.source_datenum = source_file.datenum;
    cache_info.extraction_version = extraction_version;
    cache_info.sample_count = size(v, 2);
    cache_info.created_at = datetime('now');

    save_aligned_vu_cache(aligned_data_path, v, u_log, t, cache_info);
    cache_status = sprintf('rebuilt %s', aligned_data_path);
end

function load_path = find_existing_aligned_data(aligned_data_path, legacy_aligned_data_path)
    load_path = '';
    if exist(aligned_data_path, 'file') == 2
        load_path = aligned_data_path;
    elseif exist(legacy_aligned_data_path, 'file') == 2
        load_path = legacy_aligned_data_path;
    end
end

function [v, u_log, t, cache_info] = load_aligned_vu_cache(cache_path)
    C = load(cache_path);

    if isfield(C, 'v')
        v = C.v;
    elseif isfield(C, 'y')
        v = C.y;
    else
        error('对齐数据缓存缺少 v/y 字段: %s', cache_path);
    end

    if ~isfield(C, 'u_log') || ~isfield(C, 't')
        error('对齐数据缓存缺少 u_log/t 字段: %s', cache_path);
    end

    u_log = C.u_log;
    t = C.t;
    if isfield(C, 'cache_info')
        cache_info = C.cache_info;
    else
        cache_info = struct();
    end
end

function save_aligned_vu_cache(cache_path, v, u_log, t, cache_info)
    cache_dir = fileparts(cache_path);
    if exist(cache_dir, 'dir') ~= 7
        mkdir(cache_dir);
    end

    y = v; %#ok<NASGU> 兼容旧 plot/cache 读取代码。
    save(cache_path, 'v', 'y', 'u_log', 't', 'cache_info');
end

function out = run_dp_series(v, B, umin, umax, use_restoring, opts)
    % 对每个采样点运行 direction-preserving LPCA。
    % 输出 u 的列数和 v 一致；errout/lambda 直接记录求解器返回值。
    n = size(v, 2);
    u = nan(numel(umin), n);
    errout = zeros(1, n);
    lambda = nan(1, n);

    for k = 1:n
        v_sample = v(:, k);
        try
            [~, uk, errout(k), lambda(k)] = evalc( ...
                'DP_LPCA(v_sample, B, umin, umax, 100, opts);');

            uk = clamp_cols(uk(:), umin, umax);

            if use_restoring
                uk = restoring_cpp(B, uk, umin, umax);
            end

            u(:, k) = uk;
        catch
            errout(k) = -999;
        end
    end

    out = make_result(u, errout, lambda);
end

function out = run_dpscaled_series(v, B, umin, umax, use_restoring, opts)
    % 对每个采样点运行 DPscaled_LPCA。
    % rho 是 DPscaled 内部沿目标方向的投影比例。
    n = size(v, 2);
    u = nan(numel(umin), n);
    errout = zeros(1, n);
    rho = nan(1, n);

    for k = 1:n
        v_sample = v(:, k);
        try
            [~, uk, ~, errout(k), rho(k)] = evalc( ...
                'DPscaled_LPCA(v_sample, B, umin, umax, 100, opts);');

            uk = clamp_cols(uk(:), umin, umax);

            if use_restoring
                uk = restoring_cpp(B, uk, umin, umax);
            end

            u(:, k) = uk;
        catch
            errout(k) = -999;
        end
    end

    out = make_result(u, errout, rho);
end

function u_full = add_servo_solution(u_force, servo_cols, u_servo)
    u_full = u_force;
    u_full(servo_cols, :) = u_full(servo_cols, :) + u_servo;
end

function print_summary(results, B_eval, v_eval, u_reference, umin, umax, u_reference_label)
    names = fieldnames(results);
    fprintf('\n执行器输出参考: %s\n', u_reference_label);
    fprintf('\n%-22s %10s %10s %12s %12s %10s %10s\n', ...
        'method', 'rms_u_ref', 'max_u_ref', 'rms_Bu_v', 'max_Bu_v', 'errout', 'sat');

    for i = 1:numel(names)
        name = names{i};
        u = results.(name).u;
        good = all(isfinite(u), 1);
        du = u(:, good) - u_reference(:, good);
        dv = B_eval * u(:, good) - v_eval(:, good);
        sat = any(abs(u - umin) < 1e-6 | abs(u - umax) < 1e-6, 1) & good;

        fprintf('%-22s %10.4g %10.4g %12.4g %12.4g %10d %10d\n', ...
            name, rms_vec(du), max(abs(du), [], 'all'), ...
            rms_vec(dv), max(vecnorm(dv, 2, 1)), ...
            nnz(results.(name).errout ~= 0), nnz(sat));
    end
end

function print_difference_summary(results, reference_name, B_eval, v_eval, umin, umax)
    % B_norm/B_par 坐标下的算法差异表。
    % 参考量用 inv_Bpar，因为它是 B_norm 删除零行后的直接伪逆解。
    names = fieldnames(results);
    u_ref = results.(reference_name).u;

    fprintf('\nB_norm allocation differences, reference = %s\n', reference_name);
    fprintf('%-34s %12s %12s %12s %12s %10s %10s\n', ...
        'method', 'rms_du_ref', 'max_du_ref', 'rms_Bu_v', 'max_Bu_v', 'errout', 'sat');

    for i = 1:numel(names)
        name = names{i};
        u = results.(name).u;
        good = all(isfinite(u), 1) & all(isfinite(u_ref), 1);

        du_ref = u(:, good) - u_ref(:, good);
        dv = B_eval * u(:, good) - v_eval(:, good);
        sat = any(abs(u - umin) < 1e-6 | abs(u - umax) < 1e-6, 1) & good;

        fprintf('%-34s %12.4g %12.4g %12.4g %12.4g %10d %10d\n', ...
            name, rms_vec(du_ref), max(abs(du_ref), [], 'all'), ...
            rms_vec(dv), max(vecnorm(dv, 2, 1)), ...
            nnz(results.(name).errout ~= 0), nnz(sat));
    end
end

function [v, u_log, t] = load_log_vu(mat_path)
    % 从 ulog 解析好的 MAT 里只取控制分配需要的量：
    %   v     = [vehicle_torque_setpoint.xyz; vehicle_thrust_setpoint.xyz]
    %   u_log = [actuator_motors.control0; actuator_servos.control0..5]
    S = load(mat_path);
    topics = S.log.data;

    torque_tbl = topics.vehicle_torque_setpoint_0;
    thrust_tbl = topics.vehicle_thrust_setpoint_0;
    motors_tbl = topics.actuator_motors_0;
    servos_tbl = topics.actuator_servos_0;

    t_motors = table_var(motors_tbl, {'timestamp'}) * 1e-6;
    t_control = t_motors;
    if any(strcmp(motors_tbl.Properties.VariableNames, 'timestamp_sample'))
        ts = table_var(motors_tbl, {'timestamp_sample'}) * 1e-6;
        valid = isfinite(ts) & ts > 0;
        t_control(valid) = ts(valid);
    end

    motor0 = indexed_vars(motors_tbl, 'control', 1);
    servos = indexed_vars(servos_tbl, 'control', 6);
    servos = align_to( ...
        table_var(servos_tbl, {'timestamp'}) * 1e-6, ...
        servos, t_motors, 'nearest', 0.03);

    torque = align_to( ...
        table_var(torque_tbl, {'timestamp'}) * 1e-6, ...
        xyz_vars(torque_tbl), t_control, 'nearest', 0.03);
    thrust = align_to( ...
        table_var(thrust_tbl, {'timestamp'}) * 1e-6, ...
        xyz_vars(thrust_tbl), t_control, 'previous', 0.03);

    Y = [torque, thrust];
    U = [motor0, servos];
    ok = all(isfinite(Y), 2) & all(isfinite(U), 2);

    v = Y(ok, :)';
    u_log = U(ok, :)';
    t = t_motors(ok)' - t_motors(find(ok, 1, 'first'));
end

function values_q = align_to(t, values, tq, method, max_dt_s)
    t = t(:);
    tq = tq(:);
    values = double(values);
    valid = isfinite(t) & all(isfinite(values), 2);
    [t, order] = sort(t(valid));
    values = values(valid, :);
    values = values(order, :);
    [t, unique_idx] = unique(t, 'stable');
    values = values(unique_idx, :);

    idx = round(interp1(t, (1:numel(t))', tq, method, NaN));
    values_q = nan(numel(tq), size(values, 2));
    ok = isfinite(idx) & idx >= 1 & idx <= numel(t);
    values_q(ok, :) = values(idx(ok), :);

    dt = nan(numel(tq), 1);
    dt(ok) = tq(ok) - t(idx(ok));
    if strcmp(method, 'nearest')
        dt = abs(dt);
    end
    values_q(abs(dt) > max_dt_s, :) = NaN;
end

function xyz = xyz_vars(tbl)
    xyz = [
        table_var(tbl, {'xyz_0_', 'xyz_0', 'xyz[0]', 'x'}), ...
        table_var(tbl, {'xyz_1_', 'xyz_1', 'xyz[1]', 'y'}), ...
        table_var(tbl, {'xyz_2_', 'xyz_2', 'xyz[2]', 'z'})
    ];
end

function values = indexed_vars(tbl, prefix, count)
    values = nan(height(tbl), count);
    for i = 1:count
        idx = i - 1;
        values(:, i) = table_var(tbl, {
            sprintf('%s_%d_', prefix, idx), ...
            sprintf('%s_%d', prefix, idx), ...
            sprintf('%s[%d]', prefix, idx)
        });
    end
end

function value = table_var(tbl, candidates)
    vars = tbl.Properties.VariableNames;
    desc = tbl.Properties.VariableDescriptions;
    for i = 1:numel(candidates)
        idx = find(strcmp(vars, candidates{i}), 1);
        if isempty(idx)
            idx = find(strcmp(desc, candidates{i}), 1);
        end
        if ~isempty(idx)
            value = double(tbl{:, idx});
            value = value(:);
            return;
        end
    end
    error('表中没有字段: %s', strjoin(candidates, ', '));
end

function mix = px4_geninv(G)
    G = single(G);
    [m, n] = size(G);

    if m <= n
        A = G * G';
        [L, r] = full_rank_cholesky_px4(A);
        X = rank_inverse_px4(L' * L, r);
        mix = G' * (L * (X * X * L'));
    else
        A = G' * G;
        [L, r] = full_rank_cholesky_px4(A);
        X = rank_inverse_px4(L' * L, r);
        mix = (L * (X * X * L')) * G';
    end

    mix = double(mix);
end

function scale = px4_shc09_scale_from_mix(mix)
    scale = ones(6, 1);
    fz_values = abs(mix(:, 6));
    active = fz_values > eps('single');
    scale(6) = sum(fz_values(active)) / nnz(active);
    scale(4) = scale(6);
    scale(5) = scale(6);
end

function [L, r] = full_rank_cholesky_px4(A)
    n = size(A, 1);
    tol = single(n) * eps('single') * max(diag(A));
    L = zeros(n, n, 'single');
    r = 0;

    for k = 1:n
        col = r + 1;
        if r == 0
            L(k:n, col) = A(k:n, k);
        else
            L(k:n, col) = A(k:n, k) - L(k:n, 1:r) * L(k, 1:r)';
        end
        if L(k, col) > tol
            L(k, col) = sqrt(L(k, col));
            if k < n
                L(k+1:n, col) = L(k+1:n, col) / L(k, col);
            end
            r = r + 1;
        end
    end
end

function X = rank_inverse_px4(A, r)
    X = zeros(size(A), 'single');
    if r > 0
        X(1:r, 1:r) = inv(A(1:r, 1:r));
    end
end

function u = clamp_cols(u, umin, umax)
    u = min(max(u, umin), umax);
end

function y = rms_vec(x)
    x = x(isfinite(x));
    y = sqrt(mean(x(:).^2));
end

function s = join_num(x)
    s = strjoin(string(x(:)'), ' ');
end
