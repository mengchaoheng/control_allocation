clear all;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpath(script_dir));

%% SHW09_vtol physical-B allocation model comparison
% This script tests the same 3x8 SHW09_vtol model used by the C: mixer:
%
%   MC/transition phase:
%     - bottom surfaces 1..6 use USER_OMEGA2F_MC
%     - elevons 7..8 are disabled by B(:,7:8)=0 and umin/umax=0
%
%   FW phase:
%     - bottom surfaces 1..6 use USER_OMEGA2F_FW
%     - elevons 7..8 use USER_ELEV_2_F and full physical limits
%
% It runs three cases against the same reference input y:
%   1. MC full-8: keep 8 dimensions and lock elevons 7..8 to zero
%   2. MC reduced-6: mathematically equivalent model with columns 7..8 removed
%   3. MC zeroB-8 free78: B(:,7:8)=0 but u7/u8 limits remain free
%   4. FW full-8: fixed-wing B/limits with elevons 7..8 enabled

%% User knobs
methods_to_run = {'inv', 'dir', 'dpscaled', 'prio', 'wls'};
plot_methods = methods_to_run;
enable_reachable_set_view = true;
enable_u_comparison_view = true;
enable_mc_reduced_equivalence_check = true;

% Use [] for all input. A shorter window makes iteration faster.
requested_time_window_s = [24 26];

USER_OMEGA2F_MC = 4.2;
USER_OMEGA2F_FW = 7.2;
USER_ELEV_2_F = 4.0;

%% Load reference input y
load('input.mat'); % provides v, len_command_px4, controls_delta_t_s, ...
if isfile('input.csv')
    v = readmatrix('input.csv')';
end

if size(v, 1) < 3
    error('test_shw09_vtol_modes:InputSize', 'Expected at least 3 rows in v for [roll; pitch; yaw] y input.');
end

dt = mean(controls_delta_t_s);
t_full = 0:dt:dt*(len_command_px4-1);

if isempty(requested_time_window_s)
    sample_indices = 1:len_command_px4;
elseif t_full(end) >= requested_time_window_s(1)
    sample_indices = find(t_full >= requested_time_window_s(1) & t_full <= min(requested_time_window_s(2), t_full(end)));
else
    warning('Requested window starts after input ends; using full input.');
    sample_indices = 1:len_command_px4;
end

if isempty(sample_indices)
    error('test_shw09_vtol_modes:EmptyWindow', 'No input samples selected.');
end

t = t_full(sample_indices);
y_ref = v(1:3, sample_indices);
N = size(y_ref, 2);

fprintf('Selected input window: %.3f to %.3f s, %d samples\n', t(1), t(end), N);

%% Build test cases
cases = {
    make_shw09_vtol_case('mc-full8', USER_OMEGA2F_MC, USER_OMEGA2F_FW, USER_ELEV_2_F), ...
    make_shw09_vtol_case('mc-reduced6', USER_OMEGA2F_MC, USER_OMEGA2F_FW, USER_ELEV_2_F), ...
    make_shw09_vtol_case('mc-zeroB-free78', USER_OMEGA2F_MC, USER_OMEGA2F_FW, USER_ELEV_2_F), ...
    make_shw09_vtol_case('fw-full8', USER_OMEGA2F_MC, USER_OMEGA2F_FW, USER_ELEV_2_F)
};

tie_opts.tie_rel_tol = 1e-5;
tie_opts.tie_abs_tol = 1e-6;
tie_opts.zero_tie_abs_tol = 3e-5;

%% Print model summary
fprintf('\n=== SHW09_vtol MC model ===\n');
print_model(cases{1}.B, cases{1}.umin, cases{1}.umax);
fprintf('\n=== SHW09_vtol MC reduced model ===\n');
print_model(cases{2}.B, cases{2}.umin, cases{2}.umax);
fprintf('\n=== SHW09_vtol MC zero-B free-78 model ===\n');
print_model(cases{3}.B, cases{3}.umin, cases{3}.umax);
fprintf('\n=== SHW09_vtol FW model ===\n');
print_model(cases{4}.B, cases{4}.umin, cases{4}.umax);
print_reachable_summary(cases{1}.B, cases{1}.umin, cases{1}.umax, 'MC full-8');
print_reachable_summary(cases{2}.B, cases{2}.umin, cases{2}.umax, 'MC reduced-6');
print_reachable_summary(cases{3}.B, cases{3}.umin, cases{3}.umax, 'MC zeroB-free78');
print_reachable_summary(cases{4}.B, cases{4}.umin, cases{4}.umax, 'FW full-8');

%% Simulate
results = cell(size(cases));
tic;
for case_idx = 1:numel(cases)
    results{case_idx} = simulate_case(cases{case_idx}, y_ref, t, methods_to_run, tie_opts);
end
elapsed_time = toc;

fprintf('\nMATLAB allocation test time: %.2f s\n', elapsed_time);

%% Reports
for case_idx = 1:numel(results)
    report_result(results{case_idx}, methods_to_run, t);
    report_method_differences(results{case_idx}, methods_to_run);
end

if enable_mc_reduced_equivalence_check
    compare_mc_full8_vs_reduced6(results{1}, results{2}, methods_to_run, t);
    compare_mc_locked_vs_free78(results{1}, results{3}, methods_to_run, t);
end

%% Plots and save
if enable_reachable_set_view
    plot_qcat_reachable_views(results);
end

if enable_u_comparison_view
    for plot_idx = 1:numel(plot_methods)
        plot_u_comparison_by_model(results, plot_methods{plot_idx}, t);
    end
end

save('test_shw09_vtol_modes_results.mat', 'results', 'cases', 'methods_to_run', ...
     'plot_methods', 'y_ref', 't', 'tie_opts', ...
     'USER_OMEGA2F_MC', 'USER_OMEGA2F_FW', 'USER_ELEV_2_F', ...
     'enable_reachable_set_view', 'enable_u_comparison_view', 'enable_mc_reduced_equivalence_check');

fprintf('\nSaved test_shw09_vtol_modes_results.mat\n');

function aircraft = make_shw09_vtol_case(mode, omega2f_mc, omega2f_fw, elevon_2_f)
    mc = make_shw09_vtol_model('mc', omega2f_mc, omega2f_fw, elevon_2_f);
    fw = make_shw09_vtol_model('fw', omega2f_mc, omega2f_fw, elevon_2_f);

    aircraft.name = ['SHW09-vtol-' mode];
    aircraft.mode = mode;
    aircraft.phase = string(mode);
    aircraft.plot_num_outputs = 8;

    if strcmp(mode, 'mc-full8')
        aircraft.B = mc.B;
        aircraft.umin = mc.umin;
        aircraft.umax = mc.umax;
        aircraft.plot_active_idx = 1:8;

    elseif strcmp(mode, 'mc-reduced6')
        [B_reduced, umin_reduced, umax_reduced, active_idx] = reachable_view_model(mc.B, mc.umin, mc.umax);
        aircraft.B = B_reduced;
        aircraft.umin = umin_reduced;
        aircraft.umax = umax_reduced;
        aircraft.plot_active_idx = active_idx;

    elseif strcmp(mode, 'mc-zeroB-free78')
        mc_free = make_shw09_vtol_model('mc-free78', omega2f_mc, omega2f_fw, elevon_2_f);
        aircraft.B = mc_free.B;
        aircraft.umin = mc_free.umin;
        aircraft.umax = mc_free.umax;
        aircraft.plot_active_idx = 1:8;

    elseif strcmp(mode, 'fw-full8')
        aircraft.B = fw.B;
        aircraft.umin = fw.umin;
        aircraft.umax = fw.umax;
        aircraft.plot_active_idx = 1:8;

    else
        error('Unknown SHW09_vtol test mode: %s', mode);
    end
end

function model = make_shw09_vtol_model(phase, omega2f_mc, omega2f_fw, elevon_2_f)
    I_x = 0.050636;
    I_y = 0.042954;
    I_z = 0.012668;
    L_1 = 0.42;
    L_2 = 0.073699;
    d = 60*pi/180;

    elevons_enabled = strcmp(phase, 'fw');
    elevons_locked = strcmp(phase, 'mc');

    if elevons_enabled
        k = omega2f_fw;
    else
        k = omega2f_mc;
    end

    B = zeros(3, 8);
    B(1,1:6) = [-L_1, -cos(d)*L_1, cos(d)*L_1, L_1, cos(d)*L_1, -cos(d)*L_1] * k / I_x;
    B(2,1:6) = [0, sin(d)*L_1, sin(d)*L_1, 0, -sin(d)*L_1, -sin(d)*L_1] * k / I_y;
    B(3,1:6) = ones(1,6) * L_2 * k / I_z;

    if elevons_enabled
        B(3,7) =  0.5 * L_2 * elevon_2_f / I_z;
        B(3,8) = -0.5 * L_2 * elevon_2_f / I_z;
    end

    ulim = 0.6981;
    umin = ones(8,1) * -ulim;
    umax = ones(8,1) *  ulim;

    if elevons_locked
        umin(7:8) = 0;
        umax(7:8) = 0;
    end

    model.B = B;
    model.umin = umin;
    model.umax = umax;
end

function result = simulate_case(aircraft, y_ref, t, methods_to_run, tie_opts)
    [k, N] = size(y_ref);
    m = size(aircraft.B, 2);

    result = aircraft;
    result.t = t;
    result.y_ref = y_ref;
    result.B_hist = zeros(k, m, N);
    result.umin_hist = zeros(m, N);
    result.umax_hist = zeros(m, N);
    result.phase = strings(1, N);

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};
        result.alloc.(method).u_raw = nan(m, N);
        result.alloc.(method).u = nan(m, N);
        result.alloc.(method).y_achieved = nan(k, N);
        result.alloc.(method).ok = false(1, N);
    end

    for idx = 1:N
        [B, umin, umax, phase] = model_at_time(aircraft);
        y = y_ref(:, idx);

        result.B_hist(:,:,idx) = B;
        result.umin_hist(:,idx) = umin;
        result.umax_hist(:,idx) = umax;
        result.phase(idx) = phase;

        for method_idx = 1:numel(methods_to_run)
            method = methods_to_run{method_idx};
            [u_raw, u, ok] = run_allocator(method, y, B, umin, umax, tie_opts);

            result.alloc.(method).u_raw(:,idx) = u_raw;
            result.alloc.(method).u(:,idx) = u;
            result.alloc.(method).ok(idx) = ok;

            if ok
                result.alloc.(method).y_achieved(:,idx) = B * u;
            end
        end
    end
end

function [B, umin, umax, phase] = model_at_time(aircraft)
    B = aircraft.B;
    umin = aircraft.umin;
    umax = aircraft.umax;
    phase = aircraft.phase;
end

function [u_raw, u, ok] = run_allocator(method, y, B, umin, umax, tie_opts)
    m = size(B, 2);
    k = size(B, 1);
    u_raw = nan(m, 1);
    u = nan(m, 1);
    ok = true;

    try
        switch method
            case 'inv'
                u_raw = pinv(B) * y;
                u = clamp_u(u_raw, umin, umax);

            case 'dir'
                [u_tmp, ~, ~] = DP_LPCA(y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = restoring_cpp(B, u_raw, umin, umax);

            case 'dpscaled'
                [u_tmp, ~, ~, ~] = DPscaled_LPCA(y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = restoring_cpp(B, u_raw, umin, umax);

            case 'prio'
                higher = zeros(k, 1);
                [u_tmp, ~, ~] = DP_LPCA_prio(higher, y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = restoring_cpp(B, u_raw, umin, umax);

            case 'wls'
                Wv = eye(k);
                Wu = eye(m);
                ud = zeros(m,1);
                gam = 1e6;
                W0 = zeros(m,1);
                u_tmp = wls_alloc_gen(B, y, umin, umax, Wv, Wu, ud, gam, zeros(m,1), W0, 100, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = u_raw;

            otherwise
                error('Unknown method: %s', method);
        end
    catch
        ok = false;
    end
end

function u = clamp_u(u, umin, umax)
    u = min(max(u(:), umin), umax);
end

function report_result(result, methods_to_run, t)
    fprintf('\n=== %s ===\n', result.name);
    phase_change = find(result.phase(2:end) ~= result.phase(1:end-1), 1, 'first') + 1;

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};
        data = result.alloc.(method);
        ok = data.ok;

        if ~any(ok)
            fprintf('[%s] no valid samples\n', method);
            continue;
        end

        err = data.y_achieved(:,ok) - result.y_ref(:,ok);
        err_norm = vecnorm(err, 2, 1);
        u_abs_max = max(abs(data.u(:,ok)), [], 'all');
        saturation_fraction = mean(any(abs(data.u(:,ok) - result.umin_hist(:,ok)) < 1e-6 | ...
                                       abs(data.u(:,ok) - result.umax_hist(:,ok)) < 1e-6, 1));

        err_rms = sqrt(mean(err_norm.^2));

        fprintf('[%s] ok=%d/%d, max|err|=%.6g, rms|err|=%.6g, max|u|=%.6g, sat=%.2f%%', ...
                method, nnz(ok), numel(ok), max(err_norm), err_rms, u_abs_max, 100*saturation_fraction);

        if ~isempty(phase_change) && phase_change > 1 && phase_change <= numel(t)
            du = data.u(:,phase_change) - data.u(:,phase_change-1);
            fprintf(', switch du_inf=%.6g at %.3f s', max(abs(du)), t(phase_change));
        end

        fprintf('\n');
    end
end

function report_method_differences(result, methods_to_run)
    fprintf('\n--- method differences in %s ---\n', result.name);

    baseline = 'inv';
    if ~isfield(result.alloc, baseline)
        baseline = methods_to_run{1};
    end

    base_data = result.alloc.(baseline);

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};

        if strcmp(method, baseline) || ~isfield(result.alloc, method)
            continue;
        end

        data = result.alloc.(method);
        ok = base_data.ok & data.ok;

        if ~any(ok)
            fprintf('[%s vs %s] no common valid samples\n', method, baseline);
            continue;
        end

        du = max(abs(data.u(:, ok) - base_data.u(:, ok)), [], 'all');
        dy = max(abs(data.y_achieved(:, ok) - base_data.y_achieved(:, ok)), [], 'all');
        fprintf('[%s vs %s] max|du|=%.6g, max|d(Bu)|=%.6g\n', method, baseline, du, dy);
    end

    for i = 1:numel(methods_to_run)
        for j = i+1:numel(methods_to_run)
            method_i = methods_to_run{i};
            method_j = methods_to_run{j};

            if ~isfield(result.alloc, method_i) || ~isfield(result.alloc, method_j)
                continue;
            end

            data_i = result.alloc.(method_i);
            data_j = result.alloc.(method_j);
            ok = data_i.ok & data_j.ok;

            if ~any(ok)
                continue;
            end

            du = max(abs(data_i.u(:, ok) - data_j.u(:, ok)), [], 'all');

            if du > 1e-8
                dy = max(abs(data_i.y_achieved(:, ok) - data_j.y_achieved(:, ok)), [], 'all');
                fprintf('[%s vs %s] pairwise max|du|=%.6g, max|d(Bu)|=%.6g\n', ...
                        method_i, method_j, du, dy);
            end
        end
    end
end

function compare_mc_full8_vs_reduced6(full_result, reduced_result, methods_to_run, t)
    active_idx = reduced_result.plot_active_idx;
    inactive_idx = setdiff(1:full_result.plot_num_outputs, active_idx);
    N = size(full_result.y_ref, 2);
    m_full = full_result.plot_num_outputs;

    fprintf('\n=== MC full-8 vs reduced-%d equivalence ===\n', numel(active_idx));
    fprintf('active u = ['); fprintf(' %d', active_idx); fprintf(' ], inactive u = ['); fprintf(' %d', inactive_idx); fprintf(' ]\n');

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};

        if ~isfield(full_result.alloc, method) || ~isfield(reduced_result.alloc, method)
            continue;
        end

        data_full = full_result.alloc.(method);
        data_reduced = reduced_result.alloc.(method);
        u_reduced_as_full = embed_u_for_plot(reduced_result, method, m_full);

        ok = data_full.ok & data_reduced.ok;

        if ~any(ok)
            fprintf('[%s] no common valid samples\n', method);
            continue;
        end

        du_all = max(abs(data_full.u(:, ok) - u_reduced_as_full(:, ok)), [], 'all');
        du_active = max(abs(data_full.u(active_idx, ok) - u_reduced_as_full(active_idx, ok)), [], 'all');
        dy = max(abs(data_full.y_achieved(:, ok) - data_reduced.y_achieved(:, ok)), [], 'all');

        if isempty(inactive_idx)
            u_inactive_full = 0;
        else
            u_inactive_full = max(abs(data_full.u(inactive_idx, ok)), [], 'all');
        end

        [~, worst_idx] = max(vecnorm(data_full.u(:, ok) - u_reduced_as_full(:, ok), Inf, 1));
        ok_indices = find(ok);
        worst_sample = ok_indices(worst_idx);

        fprintf('[%s] samples=%d, max|u8-u6_embed|=%.6g, active=%.6g, max|dy|=%.6g, max inactive full=%.6g, worst t=%.3f s\n', ...
                method, nnz(ok), du_all, du_active, dy, u_inactive_full, t(worst_sample));
    end
end

function compare_mc_locked_vs_free78(locked_result, free_result, methods_to_run, t)
    inactive_idx = 7:8;
    active_idx = 1:6;

    fprintf('\n=== MC locked-78 vs zeroB-free78 comparison ===\n');
    fprintf('B(:,7:8)=0 in both cases; locked has u7/u8 limits 0..0, free keeps original limits.\n');

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};

        if ~isfield(locked_result.alloc, method) || ~isfield(free_result.alloc, method)
            continue;
        end

        locked = locked_result.alloc.(method);
        free = free_result.alloc.(method);
        ok = locked.ok & free.ok;

        if ~any(ok)
            fprintf('[%s] no common valid samples\n', method);
            continue;
        end

        du_all = max(abs(locked.u(:, ok) - free.u(:, ok)), [], 'all');
        du_active = max(abs(locked.u(active_idx, ok) - free.u(active_idx, ok)), [], 'all');
        du_free_inactive = max(abs(free.u(inactive_idx, ok)), [], 'all');
        dy = max(abs(locked.y_achieved(:, ok) - free.y_achieved(:, ok)), [], 'all');

        [~, worst_idx] = max(vecnorm(free.u(inactive_idx, ok), Inf, 1));
        ok_indices = find(ok);
        worst_sample = ok_indices(worst_idx);

        fprintf('[%s] samples=%d, max|u_locked-u_free|=%.6g, active=%.6g, max|u7/u8 free|=%.6g, max|dy|=%.6g, worst t=%.3f s\n', ...
                method, nnz(ok), du_all, du_active, du_free_inactive, dy, t(worst_sample));
    end
end

function print_model(B, umin, umax)
    fprintf('B =\n');
    disp(B);
    fprintf('umin = ['); fprintf(' %.6f', umin); fprintf(' ]^T\n');
    fprintf('umax = ['); fprintf(' %.6f', umax); fprintf(' ]^T\n');
end

function vertices = reachable_vertices(B, umin, umax)
    m = length(umin);
    vertices = zeros(size(B,1), 2^m);

    for i = 0:(2^m - 1)
        u = umin;

        for j = 1:m
            if bitget(i, j)
                u(j) = umax(j);
            end
        end

        vertices(:, i+1) = B * u;
    end
end

function print_reachable_summary(B, umin, umax, name)
    V = reachable_vertices(B, umin, umax);
    fprintf('\nReachable vertex summary %s:\n', name);
    fprintf('  roll  min/max: %.6g %.6g\n', min(V(1,:)), max(V(1,:)));
    fprintf('  pitch min/max: %.6g %.6g\n', min(V(2,:)), max(V(2,:)));
    fprintf('  yaw   min/max: %.6g %.6g\n', min(V(3,:)), max(V(3,:)));
end

function plot_qcat_reachable_views(results)
    % Reference style: test_qcat_vview.m uses plim=[umin umax] and P=pinv(B)
    % with QCAT/vview.
    mc_reduced = find_result_by_mode(results, 'mc-reduced6');
    fw_full = find_result_by_mode(results, 'fw-full8');
    plot_qcat_reachable_single(mc_reduced.B, mc_reduced.umin, mc_reduced.umax, 'SHW09-vtol', 'MC reduced-6 B/limits');
    plot_qcat_reachable_single(fw_full.B, fw_full.umin, fw_full.umax, 'SHW09-vtol', 'FW full-8 B/limits');
end

function ratio = plot_qcat_reachable_single(B, umin, umax, result_name, phase_name)
    [B_view, umin_view, umax_view, active_idx] = reachable_view_model(B, umin, umax);
    plim = [umin_view umax_view];
    figure('Name', ['QCAT reachable view ' result_name ' ' phase_name]);

    if exist('vview', 'file') ~= 2
        error('test_shw09_vtol_modes:MissingQCAT', 'QCAT vview was not found on the MATLAB path.');
    end

    ratio = vview(B_view, plim, pinv(B_view));
    title(sprintf('%s (%s), active u=[%s], ratio=%.4g', ...
          result_name, phase_name, sprintf('%d ', active_idx), ratio), 'Interpreter', 'none');
end

function [B_view, umin_view, umax_view, active_idx] = reachable_view_model(B, umin, umax)
    tol = 1e-9;
    active_idx = find((abs(umax - umin) > tol) & (vecnorm(B, 2, 1)' > tol));

    if isempty(active_idx)
        error('test_shw09_vtol_modes:EmptyReachableSet', 'No active actuators for reachable-set view.');
    end

    B_view = B(:, active_idx);
    umin_view = umin(active_idx);
    umax_view = umax(active_idx);
end

function plot_u_comparison_by_model(results, method, t)
    if isempty(results)
        return;
    end

    if ~isfield(results{1}.alloc, method)
        warning('plot method %s not found', method);
        return;
    end

    n_act = max(cellfun(@(r) get_plot_num_outputs(r, method), results));
    case_names = cellfun(@case_plot_label, results, 'UniformOutput', false);
    colors = [0.0000 0.4470 0.7410;
              0.8500 0.3250 0.0980;
              0.4660 0.6740 0.1880;
              0.4940 0.1840 0.5560];
    line_styles = {'-', '--', '-.', ':'};
    line_widths = [1.15, 1.05, 1.2, 1.2];

    figure('Name', ['SHW09_vtol u comparison ' method]);
    tiledlayout(n_act, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for act_idx = 1:n_act
        nexttile;
        hold on;

        for case_idx = 1:numel(results)
            result = results{case_idx};

            if ~isfield(result.alloc, method)
                continue;
            end

            data = embed_u_for_plot(result, method, n_act);
            style_idx = min(case_idx, numel(line_styles));
            width_idx = min(case_idx, numel(line_widths));
            plot(t, data(act_idx, :), line_styles{style_idx}, ...
                 'LineWidth', line_widths(width_idx), 'Color', colors(case_idx, :), ...
                 'DisplayName', case_names{case_idx});
        end

        grid on;
        ylabel(sprintf('u_%d', act_idx));

        if act_idx == 1
            title(sprintf('Same y input, different B/limits allocation result (%s)', method), 'Interpreter', 'none');
            legend('Location', 'eastoutside', 'Interpreter', 'none');
        end

        if act_idx == n_act
            xlabel('time [s]');
        end
    end
end

function label = case_plot_label(result)
    if strcmp(result.mode, 'mc-full8')
        label = 'MC full-8, u7/u8 locked';
    elseif strcmp(result.mode, 'mc-reduced6')
        label = 'MC reduced-6';
    elseif strcmp(result.mode, 'mc-zeroB-free78')
        label = 'MC zeroB-8, u7/u8 free';
    elseif strcmp(result.mode, 'fw-full8')
        label = 'FW full-8';
    else
        label = result.name;
    end
end

function result = find_result_by_mode(results, mode)
    for case_idx = 1:numel(results)
        if strcmp(results{case_idx}.mode, mode)
            result = results{case_idx};
            return;
        end
    end

    error('test_shw09_vtol_modes:MissingResultMode', 'Could not find result mode %s.', mode);
end

function n = get_plot_num_outputs(result, method)
    if isfield(result, 'plot_num_outputs')
        n = result.plot_num_outputs;
    elseif isfield(result.alloc, method)
        n = size(result.alloc.(method).u, 1);
    else
        n = 0;
    end
end

function u_plot = embed_u_for_plot(result, method, n_act)
    u = result.alloc.(method).u;

    if size(u, 1) == n_act
        u_plot = u;
        return;
    end

    if isfield(result, 'plot_active_idx')
        active_idx = result.plot_active_idx;
    else
        active_idx = 1:size(u, 1);
    end

    u_plot = zeros(n_act, size(u, 2));
    u_plot(active_idx, :) = u;
end
