clearvars;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpath(script_dir));

% 这个脚本只测试 test.m 的三维力矩输入：
%   y = [Mx; My; Mz] = B*u
% 不使用 PX4 日志、不使用 Fz motor 分块，也不做 SHC09 6x7 归一化。
time_window_s = [20 40];
use_restoring = false;
simplex_backends = {'original', 'tiebreak'};

[B, umin, umax] = make_test_shc09_torque_model();
[y, t, input_status] = load_test_input_torque('input.mat', 'input.csv');

time_ok = t >= time_window_s(1) & t <= time_window_s(2);
y = y(:, time_ok);
t = t(time_ok);
sample_count = size(y, 2);

channel_labels = {'servo0', 'servo1', 'servo2', 'servo3', 'servo4', 'servo5'};
axis_labels = {'Mx', 'My', 'Mz'};
zero_errout = zeros(1, sample_count);
unit_scale = ones(1, sample_count);

% inv 是本测试里的参考解：直接 pinv(B)*y 后限幅。
u_inv = clamp_cols(pinv(B) * y, umin, umax);
if use_restoring
    u_inv = apply_restoring_series(B, u_inv, umin, umax);
end

results = struct();
results.inv = make_result(u_inv, zero_errout, unit_scale);

for i = 1:numel(simplex_backends)
    backend = simplex_backends{i};
    simplex_call_opts = struct('simplex_backend', backend);

    dir_name = ['pca_dir_' backend];
    scaled_name = ['pca_dpscaled_' backend];

    results.(dir_name) = run_dp_series( ...
        y, B, umin, umax, use_restoring, simplex_call_opts);
    % results.(scaled_name) = run_dpscaled_series( ...
    %     y, B, umin, umax, use_restoring, simplex_call_opts);
end

fprintf('\n=== test.m torque input tie-break test ===\n');
fprintf('input: %s\n', input_status);
fprintf('model: test.m 15006_SHC09 torque-only B, size=%dx%d, rank=%d\n', ...
    size(B, 1), size(B, 2), rank(B));
fprintf('time window: %.6g to %.6g s, samples=%d\n', ...
    time_window_s(1), time_window_s(2), sample_count);
fprintf('use_restoring: %d\n', use_restoring);
fprintf('simplex backends: %s\n', strjoin(simplex_backends, ', '));

print_summary(results, 'inv', B, y, umin, umax);

u_reference = results.inv.u;
u_reference_label = 'inv';
plot_reference_name = 'inv';
plot_results = rmfield(results, 'inv');
plot_cache_path = fullfile(script_dir, 'cache', 'last_shc09_allocation_plot_data.mat');
save(plot_cache_path, ...
    't', 'u_reference', 'u_reference_label', 'plot_results', ...
    'channel_labels', 'B', 'y', 'axis_labels', 'plot_reference_name');
fprintf('plot cache: %s\n', plot_cache_path);

plot_all_results(t, u_reference, plot_results, channel_labels, B, y, axis_labels, '', u_reference_label);

function [B, umin, umax] = make_test_shc09_torque_model()
    % 与 test.m / make_aircraft_6() 一致的 3x6 SHC09 力矩矩阵。
    % 对输出 u 来说，同一行内各执行器的相对值更直接决定
    % “这个力矩轴由谁分担”；如果某一列同时作用多个轴，
    % 这一列内的相对值决定跨轴耦合和补偿。
    l1 = 0.292166;
    l2 = 0.073699;
    k_omega2force = 1.93;
    I_x = 0.0438;
    I_y = 0.0436;
    I_z = 0.005006;
    d = 60*pi/180;
    I = diag([I_x, I_y, I_z]);

    B = I \ [
        -l1, -l1*cos(d),  l1*cos(d),  l1,  l1*cos(d), -l1*cos(d);
           0,  l1*sin(d),  l1*sin(d),  0, -l1*sin(d), -l1*sin(d);
          l2,  l2,         l2,         l2,  l2,         l2
    ] * k_omega2force;

    ulim = 40*pi/180;
    umin = ones(6, 1) * -ulim;
    umax = ones(6, 1) *  ulim;
end

function [y, t, input_status] = load_test_input_torque(mat_path, csv_path)
    S = load(mat_path);
    y = S.v;
    input_status = mat_path;

    if isfile(csv_path)
        y = readmatrix(csv_path)';
        input_status = csv_path;
    end

    sample_count = min(S.len_command_px4, size(y, 2));
    y = y(:, 1:sample_count);

    dt = mean(S.controls_delta_t_s);
    t = 0:dt:dt*(sample_count - 1);
end

function out = run_dp_series(y, B, umin, umax, use_restoring, opts)
    sample_count = size(y, 2);
    u = nan(numel(umin), sample_count);
    errout = zeros(1, sample_count);
    lambda = nan(1, sample_count);

    for k = 1:sample_count
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
    sample_count = size(y, 2);
    u = nan(numel(umin), sample_count);
    errout = zeros(1, sample_count);
    rho = nan(1, sample_count);

    for k = 1:sample_count
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

function u = apply_restoring_series(B, u, umin, umax)
    for k = 1:size(u, 2)
        u(:, k) = restoring_cpp(B, u(:, k), umin, umax);
    end
end

function out = make_result(u, errout, scale)
    out = struct('u', u, 'errout', errout, 'scale', scale);
end

function print_summary(results, reference_name, B, y, umin, umax)
    names = fieldnames(results);
    u_ref = results.(reference_name).u;

    fprintf('\n%-26s %12s %12s %12s %12s %10s %10s\n', ...
        'method', 'rms_du_inv', 'max_du_inv', 'rms_Bu_y', 'max_Bu_y', 'errout', 'sat');

    for i = 1:numel(names)
        name = names{i};
        u = results.(name).u;
        good = all(isfinite(u), 1) & all(isfinite(u_ref), 1);

        du = u(:, good) - u_ref(:, good);
        dy = B * u(:, good) - y(:, good);
        sat = any(abs(u - umin) < 1e-6 | abs(u - umax) < 1e-6, 1) & good;

        fprintf('%-26s %12.4g %12.4g %12.4g %12.4g %10d %10d\n', ...
            name, rms_vec(du), max(abs(du), [], 'all'), ...
            rms_vec(dy), max(vecnorm(dy, 2, 1)), ...
            nnz(results.(name).errout ~= 0), nnz(sat));
    end
end

function u = clamp_cols(u, umin, umax)
    u = min(max(u, umin), umax);
end

function r = rms_vec(x)
    if isempty(x)
        r = NaN;
    else
        r = sqrt(mean(x(:).^2, 'omitnan'));
    end
end
