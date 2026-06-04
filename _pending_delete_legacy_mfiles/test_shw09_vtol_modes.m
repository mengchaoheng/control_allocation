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
% It runs the same reference input y through several B/limit choices:
%   1. MC zeroB+lock78: B(:,7:8)=0 and u7/u8 are locked at 0
%   2. MC reduced-6: columns 7..8 are removed
%   3. MC zeroB-only: B(:,7:8)=0 but u7/u8 remain free
%   4. MC lock-only: B(:,7:8) keep elevon effectiveness but limits lock them
%   5. MC complete-8: B(:,7:8) keep elevon effectiveness and u7/u8 remain free

%% User knobs
methods_to_run = {'inv', 'pca_dir', 'pca_dpscaled', 'pca_prio', 'wls'};
plot_methods = methods_to_run;
enable_reachable_set_view = false;
enable_u_comparison_view = true;
enable_case_equivalence_check = true;
print_full_model_tables = false;

case_modes_to_run = {
    'mc-full8', ...           % B(:,7:8)=0, u7/u8 locked: engineering MC disabled surfaces
    'mc-reduced6', ...        % remove u7/u8 columns: mathematical reduced model
    'mc-zeroB-free78', ...    % B(:,7:8)=0, u7/u8 free: zero effectiveness only
    'mc-elevB-lock78', ...    % B(:,7:8) nonzero, u7/u8 locked: lock only
    'mc-elevB-free78' ...     % B(:,7:8) nonzero, u7/u8 free: complete MC model
};

% Important: mc-elevB-lock78 is intentionally a stress/anti-pattern case.
% DP_LPCA only maximizes lambda; it has no secondary objective that keeps u
% small or away from actuator limits. When lambda=1 is reachable, simplex may
% return any equivalent LP vertex, often with several actuators on limits.
% restoring_cpp then tries to move along null(B) toward the full-B pseudo-
% inverse solution. In mc-elevB-lock78 that pseudo-inverse target wants to use
% u7/u8, but those channels are locked by umin=umax=0, so the restoring step is
% immediately blocked by the locked channels.

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
cases = cell(1, numel(case_modes_to_run));

for case_idx = 1:numel(case_modes_to_run)
    cases{case_idx} = make_shw09_vtol_case(case_modes_to_run{case_idx}, ...
                                           USER_OMEGA2F_MC, USER_OMEGA2F_FW, USER_ELEV_2_F);
end

tie_opts = struct(); % DP/PCA 默认使用新版 simplxuprevsol_tiebreak。


%% Print model summary
fprintf('\n=== SHW09_vtol case summary ===\n');

for case_idx = 1:numel(cases)
    print_case_summary(cases{case_idx});
    print_reachable_summary(cases{case_idx}.B, cases{case_idx}.umin, cases{case_idx}.umax, cases{case_idx}.label);

    if print_full_model_tables
        print_model(cases{case_idx}.B, cases{case_idx}.umin, cases{case_idx}.umax);
    end
end

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

if enable_case_equivalence_check
    compare_cases_against_reference(results, 'mc-full8', methods_to_run, t);
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
     'case_modes_to_run', 'enable_reachable_set_view', 'enable_u_comparison_view', ...
     'enable_case_equivalence_check', 'print_full_model_tables');

fprintf('\nSaved test_shw09_vtol_modes_results.mat\n');

function aircraft = make_shw09_vtol_case(mode, omega2f_mc, omega2f_fw, elevon_2_f)
    mc = make_shw09_vtol_model('mc', omega2f_mc, omega2f_fw, elevon_2_f);
    fw = make_shw09_vtol_model('fw', omega2f_mc, omega2f_fw, elevon_2_f);

    aircraft.name = ['15008_SHW09_vtol_' mode];
    aircraft.mode = mode;
    aircraft.label = shw09_case_label(mode);
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

    elseif strcmp(mode, 'mc-elevB-lock78')
        mc_elev_lock78 = make_shw09_vtol_model('mc-elevB-lock78', omega2f_mc, omega2f_fw, elevon_2_f);
        aircraft.B = mc_elev_lock78.B;
        aircraft.umin = mc_elev_lock78.umin;
        aircraft.umax = mc_elev_lock78.umax;
        aircraft.plot_active_idx = 1:8;

    elseif strcmp(mode, 'mc-elevB-free78')
        mc_elev_free78 = make_shw09_vtol_model('mc-elevB-free78', omega2f_mc, omega2f_fw, elevon_2_f);
        aircraft.B = mc_elev_free78.B;
        aircraft.umin = mc_elev_free78.umin;
        aircraft.umax = mc_elev_free78.umax;
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

    elevon_columns_enabled = false;
    locked_idx = [];
    k = omega2f_mc;

    switch phase
        case 'mc'
            locked_idx = 7:8;

        case 'mc-free78'
            locked_idx = [];

        case 'mc-elevB-lock78'
            elevon_columns_enabled = true;
            locked_idx = 7:8;

        case 'mc-elevB-free78'
            elevon_columns_enabled = true;
            locked_idx = [];

        case 'fw'
            elevon_columns_enabled = true;
            locked_idx = [];
            k = omega2f_fw;

        otherwise
            error('Unknown SHW09_vtol phase: %s', phase);
    end

    B = zeros(3, 8);
    B(1,1:6) = [-L_1, -cos(d)*L_1, cos(d)*L_1, L_1, cos(d)*L_1, -cos(d)*L_1] * k / I_x;
    B(2,1:6) = [0, sin(d)*L_1, sin(d)*L_1, 0, -sin(d)*L_1, -sin(d)*L_1] * k / I_y;
    B(3,1:6) = ones(1,6) * L_2 * k / I_z;

    if elevon_columns_enabled
        B(3,7) =  0.5 * L_2 * elevon_2_f / I_z;
        B(3,8) = -0.5 * L_2 * elevon_2_f / I_z;
    end

    ulim = 0.6981;
    umin = ones(8,1) * -ulim;
    umax = ones(8,1) *  ulim;

    if ~isempty(locked_idx)
        umin(locked_idx) = 0;
        umax(locked_idx) = 0;
    end

    model.B = B;
    model.umin = umin;
    model.umax = umax;
end

function label = shw09_case_label(mode)
    switch mode
        case 'mc-full8'
            label = 'MC zeroB+lock78 full-8';
        case 'mc-reduced6'
            label = 'MC reduced-6';
        case 'mc-zeroB-free78'
            label = 'MC zeroB-only free78';
        case 'mc-elevB-lock78'
            label = 'MC lock-only elevB lock78';
        case 'mc-elevB-free78'
            label = 'MC complete elevB free78';
        case 'fw-full8'
            label = 'FW full-8';
        otherwise
            label = mode;
    end
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

            case 'pca_dir'
                [u_tmp, ~, ~] = DP_LPCA(y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = restoring_cpp(B, u_raw, umin, umax);

            case 'pca_dpscaled'
                [u_tmp, ~, ~, ~] = DPscaled_LPCA(y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = restoring_cpp(B, u_raw, umin, umax);

            case 'pca_prio'
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
        locked_mask = abs(result.umax_hist(:,ok) - result.umin_hist(:,ok)) < 1e-9;
        limit_mask = abs(data.u(:,ok) - result.umin_hist(:,ok)) < 1e-6 | ...
                     abs(data.u(:,ok) - result.umax_hist(:,ok)) < 1e-6;
        saturation_fraction = mean(any(limit_mask & ~locked_mask, 1));

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

function compare_cases_against_reference(results, reference_mode, methods_to_run, t)
    ref = find_result_by_mode(results, reference_mode);
    n_act = max(cellfun(@(r) r.plot_num_outputs, results));
    active_idx = 1:6;
    elevon_idx = 7:8;

    fprintf('\n=== Same input case comparison, reference = %s ===\n', ref.label);
    fprintf('du is computed after embedding reduced models back to 8 actuator channels.\n');

    for case_idx = 1:numel(results)
        current = results{case_idx};

        if strcmp(current.mode, ref.mode)
            continue;
        end

        fprintf('\n--- %s vs %s ---\n', current.label, ref.label);

        for method_idx = 1:numel(methods_to_run)
            method = methods_to_run{method_idx};

            if ~isfield(ref.alloc, method) || ~isfield(current.alloc, method)
                continue;
            end

            ref_data = ref.alloc.(method);
            cur_data = current.alloc.(method);
            ok = ref_data.ok & cur_data.ok;

            if ~any(ok)
                fprintf('[%s] no common valid samples\n', method);
                continue;
            end

            u_ref = embed_u_for_plot(ref, method, n_act);
            u_cur = embed_u_for_plot(current, method, n_act);

            du_all = max(abs(u_cur(:, ok) - u_ref(:, ok)), [], 'all');
            du_active = max(abs(u_cur(active_idx, ok) - u_ref(active_idx, ok)), [], 'all');
            du_elevon = max(abs(u_cur(elevon_idx, ok) - u_ref(elevon_idx, ok)), [], 'all');
            dy = max(abs(cur_data.y_achieved(:, ok) - ref_data.y_achieved(:, ok)), [], 'all');

            [~, worst_idx] = max(vecnorm(u_cur(:, ok) - u_ref(:, ok), Inf, 1));
            ok_indices = find(ok);
            worst_sample = ok_indices(worst_idx);

            fprintf('[%s] samples=%d, max|du|=%.6g, active1-6=%.6g, elevon7-8=%.6g, max|d(Bu)|=%.6g, worst t=%.3f s\n', ...
                    method, nnz(ok), du_all, du_active, du_elevon, dy, t(worst_sample));
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

function print_case_summary(aircraft)
    locked_idx = find(abs(aircraft.umax - aircraft.umin) < 1e-9);
    free_idx = setdiff(1:numel(aircraft.umin), locked_idx);
    col_norm = vecnorm(aircraft.B, 2, 1);
    col7_norm = nan;
    col8_norm = nan;

    if numel(col_norm) >= 7
        col7_norm = col_norm(7);
    end

    if numel(col_norm) >= 8
        col8_norm = col_norm(8);
    end

    fprintf('\nCase: %s\n', aircraft.label);
    fprintf('  mode=%s, size(B)=%dx%d, rank(B)=%d\n', ...
            aircraft.mode, size(aircraft.B,1), size(aircraft.B,2), rank(aircraft.B));
    fprintf('  free u = ['); fprintf(' %d', free_idx); fprintf(' ], locked u = ['); fprintf(' %d', locked_idx); fprintf(' ]\n');
    fprintf('  ||B(:,7)||=%.6g, ||B(:,8)||=%.6g\n', col7_norm, col8_norm);
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
    mc_complete = find_result_by_mode(results, 'mc-elevB-free78');
    plot_qcat_reachable_single(mc_reduced.B, mc_reduced.umin, mc_reduced.umax, 'SHW09-vtol', 'MC reduced-6 B/limits');
    plot_qcat_reachable_single(mc_complete.B, mc_complete.umin, mc_complete.umax, 'SHW09-vtol', 'MC complete-8 B/limits');
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
    colors = lines(numel(results));
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
            style_idx = mod(case_idx - 1, numel(line_styles)) + 1;
            width_idx = mod(case_idx - 1, numel(line_widths)) + 1;
            plot(t, data(act_idx, :), line_styles{style_idx}, ...
                 'LineWidth', line_widths(width_idx), 'Color', colors(case_idx, :), ...
                 'DisplayName', case_names{case_idx});
        end

        grid on;
        ylabel(sprintf('u_%d', act_idx), 'Interpreter', 'none');

        if act_idx == 1
            title(sprintf('Same y input, different B/limits allocation result (%s)', method), 'Interpreter', 'none');
            legend('Location', 'best', 'Interpreter', 'none', 'NumColumns', 2);
        end

        if act_idx == n_act
            xlabel('time [s]', 'Interpreter', 'none');
        end
    end
end

function label = case_plot_label(result)
    if isfield(result, 'label')
        label = result.label;
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
