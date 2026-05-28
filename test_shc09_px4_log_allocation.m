clearvars;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpath(script_dir));

mat_path = '/Users/mch/Proj/PX4_ulog_plottools/data/06_47_22.mat';
aligned_cache_dir = fullfile(script_dir, 'cache');
plot_cache_path = fullfile(aligned_cache_dir, 'last_shc09_allocation_plot_data.mat');
force_rebuild_aligned_cache = false;
time_window_s = [8 32];      % 只跑这个时间窗口内的数据，例如 [10 20]。
use_restoring = true;

% simplex_backends 写什么就跑什么。日常使用 tiebreak：
%   - entering variable 用确定 tie-break，减少等价顶点跳变；
%   - leaving ratio 用内部 single-eps 级容差，避免把不该并列的边界误判成并列。
simplex_backends = {'tiebreak'};

% SHC09: y = B*u, y=[Mx My Mz Fx Fy Fz]', u=[motor0 servo0 ... servo5]'.
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

[y, u_reference, t, cache_status] = load_or_make_aligned_yu( ...
    mat_path, aligned_cache_dir, force_rebuild_aligned_cache);
u_reference_label = 'u_log';

t = t(:)';
time_ok = t >= time_window_s(1) & t <= time_window_s(2);
y = y(:, time_ok);
u_reference = u_reference(:, time_ok);
t = t(time_ok);

sample_count = size(y, 2);
zero_errout = zeros(1, sample_count);
unit_scale = ones(1, sample_count);

% PX4 mix normalization: mix_norm = B^+ * D^-1, 等价用 B_norm = D*B做分配。
mix_px4_raw = px4_geninv(B);
scale = px4_shc09_scale_from_mix(mix_px4_raw);
mix_px4_norm = mix_px4_raw ./ scale';
B_norm = diag(scale) * B;

active_rows = find(any(abs(B_norm) > 1e-10, 2));
B_par = B_norm(active_rows, :);
y_par = y(active_rows, :);
axis_par = axis_labels(active_rows);

results = struct();

% 1. PX4 日志对应方法：原始 B 取 geninv 后按 PX4 scale 归一化 mix。
% 这里用于验证日志中的 PX4 inv 输出，不作为后续 PCA 对比的 B。
u_px4_raw = mix_px4_norm * y;
u_px4 = clamp_cols(u_px4_raw, umin, umax);
results.px4_inv = make_result(u_px4, zero_errout, unit_scale);

% 2. 不拆力/力矩：完整 B_norm  
mix_inv_Bnorm = pinv(B_norm);
u_inv_Bnorm_raw = mix_inv_Bnorm * y;
u_inv_Bnorm = clamp_cols(u_inv_Bnorm_raw, umin, umax);
results.inv_Bnorm = make_result(u_inv_Bnorm, zero_errout, unit_scale);

% 3. 删除零行后的 B_par。B_par 只是删掉 Fx/Fy 两个全零行。
% 3.1 伪逆法。inv_Bnorm 和 inv_Bpar 应当几乎一致；
mix_inv_Bpar = pinv(B_par);
u_inv_Bpar_raw = mix_inv_Bpar * y_par;
u_inv_Bpar = clamp_cols(u_inv_Bpar_raw, umin, umax);
results.inv_Bpar = make_result(u_inv_Bpar, zero_errout, unit_scale);

if isempty(u_reference)
    u_reference = u_inv_Bpar;
end

% 4. split：先分配力，再补偿力执行器附带的力矩，最后分配剩余力矩。
%
% 这个分解依赖 SHC09 的结构：
%   y = B*u, y=[tau; f], u=[u_force; u_servo]
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
u_force(motor_col, :) = y_par(fz_row, :) / B_par(fz_row, motor_col);
u_force = clamp_cols(u_force, umin, umax);

% motor0 已经用于 Fz；剩余力矩只交给 6 个 servo。
B_torque_servo = B_par(torque_rows, servo_cols);
umin_servo = umin(servo_cols);
umax_servo = umax(servo_cols);

% Fz 的 motor0 输出也会带来 Mz，所以先从三维力矩期望里扣除这部分。
y_torque_from_motor = B_par(torque_rows, :) * u_force;
y_torque_left = y_par(torque_rows, :) - y_torque_from_motor;

% 4.1 split / inv：只对剩余三维力矩做 3x6 pinv 分配。
mix_torque_inv = pinv(B_torque_servo);
u_torque_inv_raw = mix_torque_inv * y_torque_left;
u_torque_inv = clamp_cols(u_torque_inv_raw, umin_servo, umax_servo);
u_split_inv = add_servo_solution(u_force, servo_cols, u_torque_inv);
results.split_inv = make_result(u_split_inv, zero_errout, unit_scale);

% 5. PCA/DP 方法。simplex_backends 列了几个后端，就跑几个后端。
for i = 1:numel(simplex_backends)
    backend = simplex_backends{i};
    simplex_call_opts = struct('simplex_backend', backend);

    % 不拆力/力矩：在同一个 4x7 B_par 上运行 DP / DPscaled。
    dp_name = ['pca_dir_' backend '_Bpar'];
    dps_name = ['pca_dpscaled_' backend '_Bpar'];
    results.(dp_name) = run_dp_series( ...
        y_par, B_par, umin, umax, use_restoring, simplex_call_opts);
    results.(dps_name) = run_dpscaled_series( ...
        y_par, B_par, umin, umax, use_restoring, simplex_call_opts);

    % split：Fz 仍直接算，三维剩余力矩用不同后端分配。
    split_dp_name = ['split_pca_dir_' backend];
    split_dps_name = ['split_pca_dpscaled_' backend];

    torque_dir_i = run_dp_series( ...
        y_torque_left, B_torque_servo, umin_servo, umax_servo, ...
        use_restoring, simplex_call_opts);
    torque_scaled_i = run_dpscaled_series( ...
        y_torque_left, B_torque_servo, umin_servo, umax_servo, ...
        use_restoring, simplex_call_opts);

    results.(split_dp_name) = make_result( ...
        add_servo_solution(u_force, servo_cols, torque_dir_i.u), ...
        torque_dir_i.errout, torque_dir_i.scale);
    results.(split_dps_name) = make_result( ...
        add_servo_solution(u_force, servo_cols, torque_scaled_i.u), ...
        torque_scaled_i.errout, torque_scaled_i.scale);
end

fprintf('\nB_norm rows kept for B_par: [%s] = [%s]\n', join_num(active_rows), strjoin(axis_par, ' '));
fprintf('aligned y/u cache: %s\n', cache_status);
fprintf('u reference for rms/black line: %s\n', u_reference_label);
fprintf('time window: %.6g to %.6g s, samples=%d\n', ...
    time_window_s(1), time_window_s(2), sample_count);
fprintf('simplex backends: %s\n', strjoin(simplex_backends, ', '));
print_summary(results, B_par, y_par, u_reference, umin, umax, u_reference_label);
print_difference_summary(results, 'inv_Bpar', B_par, y_par, umin, umax);

plot_reference_name = 'inv_Bpar';
save(plot_cache_path, ...
    't', 'u_reference', 'u_reference_label', 'results', ...
    'channel_labels', 'B_par', 'y_par', 'axis_par', 'plot_reference_name');
fprintf('plot cache: %s\n', plot_cache_path);

% test 跑完立即画图；之后也可以直接运行 plot_all_results 重画上次结果。
plot_all_results(t, u_reference, results, channel_labels, B_par, y_par, axis_par, 'inv_Bpar', u_reference_label);

function out = make_result(u, errout, scale)
    out = struct('u', u, 'errout', errout, 'scale', scale);
end

function [y, u_log, t, cache_status] = load_or_make_aligned_yu(mat_path, cache_dir, force_rebuild)
    % 控制分配只需要对齐后的 y 和 u_log。
    % 第一次从完整日志 MAT 提取；之后直接读取小缓存，避免反复解析 ulog 表。
    extraction_version = 1;

    [~, mat_name] = fileparts(mat_path);
    cache_path = fullfile(cache_dir, [mat_name '_aligned_yu.mat']);
    source_file = dir(mat_path);
    if isempty(source_file)
        error('找不到日志 MAT 文件: %s', mat_path);
    end

    cache_is_valid = false;
    if ~force_rebuild && exist(cache_path, 'file') == 2
        C = load(cache_path, 'y', 'u_log', 't', 'cache_info');

        cache_is_valid = isfield(C, 'y') ...
            && isfield(C, 'u_log') ...
            && isfield(C, 't') ...
            && isfield(C, 'cache_info') ...
            && strcmp(C.cache_info.source_path, mat_path) ...
            && C.cache_info.source_bytes == source_file.bytes ...
            && abs(C.cache_info.source_datenum - source_file.datenum) < 1e-9 ...
            && C.cache_info.extraction_version == extraction_version;
    end

    if cache_is_valid
        y = C.y;
        u_log = C.u_log;
        t = C.t;
        cache_status = sprintf('loaded %s', cache_path);
        return;
    end

    [y, u_log, t] = load_log_yu(mat_path);

    if exist(cache_dir, 'dir') ~= 7
        mkdir(cache_dir);
    end

    cache_info = struct();
    cache_info.source_path = mat_path;
    cache_info.source_bytes = source_file.bytes;
    cache_info.source_datenum = source_file.datenum;
    cache_info.extraction_version = extraction_version;
    cache_info.sample_count = size(y, 2);
    cache_info.created_at = datetime('now');

    save(cache_path, 'y', 'u_log', 't', 'cache_info');
    cache_status = sprintf('rebuilt %s', cache_path);
end

function out = run_dp_series(y, B, umin, umax, use_restoring, opts)
    % 对每个采样点运行 direction-preserving LPCA。
    % 输出 u 的列数和 y 一致；errout/lambda 直接记录求解器返回值。
    n = size(y, 2);
    u = nan(numel(umin), n);
    errout = zeros(1, n);
    lambda = nan(1, n);

    for k = 1:n
        y_sample = y(:, k);
        try
            [~, uk, errout(k), lambda(k)] = evalc( ...
                'DP_LPCA(y_sample, B, umin, umax, 100, opts);');

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

function out = run_dpscaled_series(y, B, umin, umax, use_restoring, opts)
    % 对每个采样点运行 DPscaled_LPCA。
    % rho 是 DPscaled 内部沿目标方向的投影比例。
    n = size(y, 2);
    u = nan(numel(umin), n);
    errout = zeros(1, n);
    rho = nan(1, n);

    for k = 1:n
        y_sample = y(:, k);
        try
            [~, uk, ~, errout(k), rho(k)] = evalc( ...
                'DPscaled_LPCA(y_sample, B, umin, umax, 100, opts);');

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

function print_summary(results, B_eval, y_eval, u_reference, umin, umax, u_reference_label)
    names = fieldnames(results);
    fprintf('\n执行器输出参考: %s\n', u_reference_label);
    fprintf('\n%-22s %10s %10s %12s %12s %10s %10s\n', ...
        'method', 'rms_u_ref', 'max_u_ref', 'rms_Bu_y', 'max_Bu_y', 'errout', 'sat');

    for i = 1:numel(names)
        name = names{i};
        u = results.(name).u;
        good = all(isfinite(u), 1);
        du = u(:, good) - u_reference(:, good);
        dy = B_eval * u(:, good) - y_eval(:, good);
        sat = any(abs(u - umin) < 1e-6 | abs(u - umax) < 1e-6, 1) & good;

        fprintf('%-22s %10.4g %10.4g %12.4g %12.4g %10d %10d\n', ...
            name, rms_vec(du), max(abs(du), [], 'all'), ...
            rms_vec(dy), max(vecnorm(dy, 2, 1)), ...
            nnz(results.(name).errout ~= 0), nnz(sat));
    end
end

function print_difference_summary(results, reference_name, B_eval, y_eval, umin, umax)
    % B_norm/B_par 坐标下的算法差异表。
    % 参考量用 inv_Bpar，因为它是 B_norm 删除零行后的直接伪逆解。
    names = fieldnames(results);
    u_ref = results.(reference_name).u;

    fprintf('\nB_norm allocation differences, reference = %s\n', reference_name);
    fprintf('%-34s %12s %12s %12s %12s %10s %10s\n', ...
        'method', 'rms_du_ref', 'max_du_ref', 'rms_Bu_y', 'max_Bu_y', 'errout', 'sat');

    for i = 1:numel(names)
        name = names{i};
        u = results.(name).u;
        good = all(isfinite(u), 1) & all(isfinite(u_ref), 1);

        du_ref = u(:, good) - u_ref(:, good);
        dy = B_eval * u(:, good) - y_eval(:, good);
        sat = any(abs(u - umin) < 1e-6 | abs(u - umax) < 1e-6, 1) & good;

        fprintf('%-34s %12.4g %12.4g %12.4g %12.4g %10d %10d\n', ...
            name, rms_vec(du_ref), max(abs(du_ref), [], 'all'), ...
            rms_vec(dy), max(vecnorm(dy, 2, 1)), ...
            nnz(results.(name).errout ~= 0), nnz(sat));
    end
end

function [y, u_log, t] = load_log_yu(mat_path)
    % 从 ulog 解析好的 MAT 里只取控制分配需要的量：
    %   y     = [vehicle_torque_setpoint.xyz; vehicle_thrust_setpoint.xyz]
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

    y = Y(ok, :)';
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
