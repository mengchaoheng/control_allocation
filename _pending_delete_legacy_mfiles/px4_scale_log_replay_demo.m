clearvars;
close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);

% Replay PX4 control-allocation inputs from a log and compare allocation
% results for different B matrices. Input and logged output are not aligned to
% each other for plotting; each uses its own timestamp sequence.

mat_path = '/Users/mch/Proj/PX4_ulog_plottools/data/15_58_23.mat';
time_window_s = [8 32];   % [] for all data
make_plots = true;
print_unitized_matrices = true;

axis_labels = {'Mx', 'My', 'Mz', 'Fx', 'Fy', 'Fz'};
channel_labels = {'motor0', 'servo0', 'servo1', 'servo2', 'servo3', 'servo4', 'servo5'};

% Edit/add entries here to compare different B matrices.
B_cases = struct([]);

B_cases(end+1).name = 'px4_status_B_motor_Mz_0p325';
B_cases(end).B = [
    0        -0.5010  -0.2505   0.2505   0.5010   0.2505  -0.2505
    0         0        0.4339   0.4339   0       -0.4339  -0.4339
    0.3250    0.2070   0.2070   0.2070   0.2070   0.2070   0.2070
    0         0        0        0        0        0        0
    0         0        0        0        0        0        0
   -6.5000    0        0        0        0        0        0
];

% B_cases(end+1).name = 'comparison_B_motor_Mz_0';
% B_cases(end).B = [
%     0        -0.5010  -0.2505   0.2505   0.5010   0.2505  -0.2505
%     0         0        0.4339   0.4339   0       -0.4339  -0.4339
%     0         0.2070   0.2070   0.2070   0.2070   0.2070   0.2070
%     0         0        0        0        0        0        0
%     0         0        0        0        0        0        0
%    -6.5000    0        0        0        0        0        0
% ];

umin = [0; -1; -1; -1; -1; -1; -1];
umax = [1;  1;  1;  1;  1;  1;  1];

% PX4 Custom effectiveness normally uses normalize_rpy=false. Keep the
% switches here so you can compare the exact scale behavior.
cfg = struct();
cfg.metric_allocation = false;
cfg.normalize_rpy = false;
cfg.update_normalization_scale = true;
cfg.use_px4_geninv = true;

[input_log, output_log] = load_log_input_output_no_output_alignment(mat_path, time_window_s);

fprintf('Loaded input samples: %d\n', numel(input_log.t));
fprintf('Loaded motor samples: %d\n', numel(output_log.t_motor));
fprintf('Loaded servo samples: %d\n', numel(output_log.t_servo));

results = struct([]);

for case_idx = 1:numel(B_cases)
    B = B_cases(case_idx).B;

    together = px4_allocate_series(B, umin, umax, input_log.c, cfg);

    c_force = zeros(size(input_log.c));
    c_force(4:6, :) = input_log.c(4:6, :);
    c_torque = zeros(size(input_log.c));
    c_torque(1:3, :) = input_log.c(1:3, :);

    force_first = px4_force_then_torque_series(B, umin, umax, c_force, c_torque, cfg);
    torque_only = px4_allocate_series(B, umin, umax, c_torque, cfg);

    results(case_idx).name = B_cases(case_idx).name;
    results(case_idx).B = B;
    results(case_idx).together = together;
    results(case_idx).force_first = force_first;
    results(case_idx).torque_only = torque_only;

    fprintf('\n%s\n', results(case_idx).name);
    fprintf('scale = [');
    fprintf(' %.6g', together.scale);
    fprintf(' ]\n');

    if print_unitized_matrices
        fprintf('D_scale = diag(scale):\n');
        disp(together.scale_matrix);
        fprintf('B_unit = D_scale * B:\n');
        disp(together.B_unit);
    end

    fprintf('rms unallocated together = %.6g\n', rms_vec(together.unallocated));
    fprintf('rms unallocated force_first = %.6g\n', rms_vec(force_first.unallocated));
    fprintf('rms unallocated torque_only torque rows = %.6g\n', rms_vec(torque_only.unallocated(1:3, :)));
end

if make_plots
    plot_input_commands(input_log, axis_labels);
    plot_logged_outputs(output_log, channel_labels);
    plot_allocated_controls(input_log.t, input_log.c, results, axis_labels);
    plot_allocation_results(input_log.t, results, channel_labels, output_log);
    
end

%% Local functions
function [input_log, output_log] = load_log_input_output_no_output_alignment(mat_path, time_window_s)
    S = load(mat_path);
    topics = S.log.data;

    torque_tbl = topics.vehicle_torque_setpoint_0;
    thrust_tbl = topics.vehicle_thrust_setpoint_0;
    motors_tbl = topics.actuator_motors_0;
    servos_tbl = topics.actuator_servos_0;

    t_torque_abs = table_var(torque_tbl, {'timestamp'}) * 1e-6;
    t_thrust_abs = table_var(thrust_tbl, {'timestamp'}) * 1e-6;
    t_motor_abs = table_var(motors_tbl, {'timestamp'}) * 1e-6;
    t_servo_abs = table_var(servos_tbl, {'timestamp'}) * 1e-6;

    t0 = min([first_finite(t_torque_abs), first_finite(t_thrust_abs), ...
              first_finite(t_motor_abs), first_finite(t_servo_abs)]);

    torque = xyz_vars(torque_tbl);
    thrust = xyz_vars(thrust_tbl);

    % This only combines the two input topics into one allocator input stream.
    % The logged actuator output is intentionally not aligned to this input.
    thrust_at_torque = align_previous(t_thrust_abs, thrust, t_torque_abs, 0.05);
    c = [torque, thrust_at_torque]';
    t_input = t_torque_abs(:)' - t0;

    ok_input = all(isfinite(c), 1) & isfinite(t_input);
    if ~isempty(time_window_s)
        ok_input = ok_input & t_input >= time_window_s(1) & t_input <= time_window_s(2);
    end

    input_log.t = t_input(ok_input);
    input_log.c = c(:, ok_input);

    motor0 = indexed_vars(motors_tbl, 'control', 1)';
    servos = indexed_vars(servos_tbl, 'control', 6)';

    t_motor = t_motor_abs(:)' - t0;
    t_servo = t_servo_abs(:)' - t0;

    ok_motor = all(isfinite(motor0), 1) & isfinite(t_motor);
    ok_servo = all(isfinite(servos), 1) & isfinite(t_servo);
    if ~isempty(time_window_s)
        ok_motor = ok_motor & t_motor >= time_window_s(1) & t_motor <= time_window_s(2);
        ok_servo = ok_servo & t_servo >= time_window_s(1) & t_servo <= time_window_s(2);
    end

    output_log.t_motor = t_motor(ok_motor);
    output_log.motor0 = motor0(:, ok_motor);
    output_log.t_servo = t_servo(ok_servo);
    output_log.servos = servos(:, ok_servo);
end

function result = px4_allocate_series(B, umin, umax, c, cfg)
    num_axes = 6;
    assert(size(B, 1) == num_axes);

    scale_info = px4_pinv_scale_info(B, cfg);

    u_raw = scale_info.mix_norm * c;
    u = clamp_cols(u_raw, umin, umax);
    allocated = scale_info.B_unit * u;

    result.mix_raw = scale_info.mix_raw;
    result.mix = scale_info.mix_norm;
    result.scale = scale_info.scale;
    result.scale_matrix = scale_info.scale_matrix;
    result.B_unit = scale_info.B_unit;
    result.u_raw = u_raw;
    result.u = u;
    result.allocated = allocated;
    result.raw_Bu = B * u;
    result.unallocated = c - allocated;
end

function scale_info = px4_pinv_scale_info(B, cfg)
    % PX4 pseudo-inverse normalization:
    %   mix_raw  = pinv(B)
    %   scale    = scale(mix_raw)
    %   mix_norm = mix_raw / scale_per_axis
    %
    % For allocation analysis it is often easier to move that same scaling
    % to the effectiveness side:
    %   allocated_control = (B * u) .* scale = diag(scale) * B * u
    % so B_unit = diag(scale) * B is the matrix PX4 effectively reports in
    % allocated_control coordinates.
    num_axes = 6;

    if cfg.use_px4_geninv
        mix_raw = px4_geninv(B);
    else
        mix_raw = pinv(B);
    end

    scale = ones(num_axes, 1);
    mix_norm = mix_raw;

    if ~cfg.metric_allocation
        if cfg.update_normalization_scale
            scale = px4_update_control_allocation_matrix_scale(mix_raw, cfg.normalize_rpy);
        end

        mix_norm = px4_normalize_control_allocation_matrix(mix_raw, scale);
    end

    scale_matrix = diag(scale);

    scale_info = struct();
    scale_info.scale = scale;
    scale_info.scale_matrix = scale_matrix;
    scale_info.B_unit = scale_matrix * B;
    scale_info.mix_raw = mix_raw;
    scale_info.mix_norm = mix_norm;
end

function result = px4_force_then_torque_series(B, umin, umax, c_force, c_torque, cfg)
    base = px4_allocate_series(B, umin, umax, zeros(size(c_force)), cfg);

    num_actuators = size(B, 2);
    sample_count = size(c_force, 2);
    umin = umin(:);
    umax = umax(:);

    % Work in PX4 allocated_control coordinates:
    %   allocated_control = (B * u) .* scale = diag(scale) * B * u.
    % This is the same B_norm coordinate used in test_shc09_px4_log_allocation.m.
    B_scaled = base.B_unit;
    torque_rows = 1:3;
    force_rows = 4:6;
    tol = 1e-10;

    force_active_rows = force_rows(any(abs(B_scaled(force_rows, :)) > tol, 2));
    force_cols = find(any(abs(B_scaled(force_active_rows, :)) > tol, 1));

    u_force_raw = zeros(num_actuators, sample_count);

    if ~isempty(force_active_rows) && ~isempty(force_cols)
        B_force = B_scaled(force_active_rows, force_cols);
        force_cmd = c_force(force_active_rows, :);

        if cfg.use_px4_geninv
            mix_force = px4_geninv(B_force);
        else
            mix_force = pinv(B_force);
        end

        u_force_raw(force_cols, :) = mix_force * force_cmd;
    end

    % Use the clamped force-actuator output for compensation. If motor thrust
    % saturates, the actual motor torque coupled into Mx/My/Mz must use the
    % saturated value, not the raw force solution.
    u_force = clamp_cols(u_force_raw, umin, umax);
    force_allocated = B_scaled * u_force;
    torque_from_force = B_scaled(torque_rows, :) * u_force;
    torque_left = c_torque(torque_rows, :) - torque_from_force;

    torque_cols = setdiff(1:num_actuators, force_cols, 'stable');
    torque_cols = torque_cols(any(abs(B_scaled(torque_rows, torque_cols)) > tol, 1));

    u_torque_raw = zeros(num_actuators, sample_count);

    if ~isempty(torque_cols)
        B_torque = B_scaled(torque_rows, torque_cols);

        if cfg.use_px4_geninv
            mix_torque = px4_geninv(B_torque);
        else
            mix_torque = pinv(B_torque);
        end

        u_torque_raw(torque_cols, :) = mix_torque * torque_left;
    end

    total_target = c_force + c_torque;
    u_raw = u_force + u_torque_raw;
    u = clamp_cols(u_raw, umin, umax);
    allocated = B_scaled * u;

    result.mix_raw = base.mix_raw;
    result.mix = base.mix;
    result.scale = base.scale;
    result.B_scaled = B_scaled;
    result.force_cols = force_cols;
    result.torque_cols = torque_cols;
    result.torque_from_force = torque_from_force;
    result.torque_left = torque_left;
    result.u_raw = u_raw;
    result.u = u;
    result.allocated = allocated;
    result.raw_Bu = B * u;
    result.unallocated = total_target - allocated;
    result.force_stage = struct( ...
        'u_raw', u_force_raw, ...
        'u', u_force, ...
        'allocated', force_allocated, ...
        'raw_Bu', B * u_force, ...
        'unallocated', c_force - force_allocated);
end

function scale = px4_update_control_allocation_matrix_scale(mix, normalize_rpy)
    FLT_EPS = double(single(eps('single')));
    scale = ones(6, 1);

    if normalize_rpy
        num_non_zero_roll = sum(abs(mix(:, 1)) > 1e-3);
        num_non_zero_pitch = sum(abs(mix(:, 2)) > 1e-3);

        roll_norm_scale = 1;
        if num_non_zero_roll > 0
            roll_norm_scale = sqrt(sum(mix(:, 1).^2) / (num_non_zero_roll / 2));
        end

        pitch_norm_scale = 1;
        if num_non_zero_pitch > 0
            pitch_norm_scale = sqrt(sum(mix(:, 2).^2) / (num_non_zero_pitch / 2));
        end

        scale(1) = max(roll_norm_scale, pitch_norm_scale);
        scale(2) = scale(1);
        scale(3) = max(mix(:, 3));
    else
        scale(1:3) = 1;
    end

    scale(6) = 1;

    for axis_idx = 2:-1:0
        j = 4 + axis_idx;  % Fx,Fy,Fz -> 4,5,6
        col_abs = abs(mix(:, j));
        num_non_zero = sum(col_abs > FLT_EPS);

        if num_non_zero > 0
            scale(j) = sum(col_abs) / num_non_zero;
        else
            scale(j) = scale(6);
        end
    end
end

function mix_norm = px4_normalize_control_allocation_matrix(mix, scale)
    FLT_EPS = double(single(eps('single')));
    mix_norm = mix;

    if scale(1) > FLT_EPS
        mix_norm(:, 1) = mix_norm(:, 1) / scale(1);
        mix_norm(:, 2) = mix_norm(:, 2) / scale(2);
    end

    if scale(3) > FLT_EPS
        mix_norm(:, 3) = mix_norm(:, 3) / scale(3);
    end

    if scale(4) > FLT_EPS
        mix_norm(:, 4) = mix_norm(:, 4) / scale(4);
        mix_norm(:, 5) = mix_norm(:, 5) / scale(5);
        mix_norm(:, 6) = mix_norm(:, 6) / scale(6);
    end

    mix_norm(abs(mix_norm) < 1e-3) = 0;
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

function values_q = align_previous(t, values, tq, max_dt_s)
    t = t(:);
    tq = tq(:);
    values = double(values);
    valid = isfinite(t) & all(isfinite(values), 2);
    [t, order] = sort(t(valid));
    values = values(valid, :);
    values = values(order, :);
    [t, unique_idx] = unique(t, 'stable');
    values = values(unique_idx, :);

    idx = round(interp1(t, (1:numel(t))', tq, 'previous', NaN));
    values_q = nan(numel(tq), size(values, 2));
    ok = isfinite(idx) & idx >= 1 & idx <= numel(t);
    values_q(ok, :) = values(idx(ok), :);

    dt = nan(numel(tq), 1);
    dt(ok) = tq(ok) - t(idx(ok));
    values_q(dt < 0 | dt > max_dt_s, :) = NaN;
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

    error('Table field not found: %s', strjoin(candidates, ', '));
end

function y = clamp_cols(x, lower, upper)
    y = min(max(x, lower), upper);
end

function t = first_finite(x)
    idx = find(isfinite(x), 1, 'first');
    if isempty(idx)
        error('timestamp vector has no finite values');
    end
    t = x(idx);
end

function y = rms_vec(x)
    x = x(isfinite(x));
    y = sqrt(mean(x(:).^2));
end

function plot_input_commands(input_log, axis_labels)
    figure('Name', 'PX4 allocation input setpoints');
    tiledlayout(6, 1, 'TileSpacing', 'compact');
    for i = 1:6
        nexttile;
        plot(input_log.t, input_log.c(i, :), 'LineWidth', 1.0);
        grid on;
        ylabel(axis_labels{i});
        if i == 1
            title('Input control_sp = [torque; thrust], own input timestamp');
        end
    end
    xlabel('input time [s]');
end

function plot_logged_outputs(output_log, channel_labels)
    figure('Name', 'PX4 logged actuator allocation output');
    tiledlayout(7, 1, 'TileSpacing', 'compact');

    nexttile;
    plot(output_log.t_motor, output_log.motor0, 'k', 'LineWidth', 1.0);
    grid on;
    ylabel(channel_labels{1});
    title('Logged actuator output, own output timestamps');

    for i = 1:6
        nexttile;
        plot(output_log.t_servo, output_log.servos(i, :), 'k', 'LineWidth', 1.0);
        grid on;
        ylabel(channel_labels{i + 1});
    end
    xlabel('output time [s]');
end

function plot_allocation_results(t, results, channel_labels, output_log)
    figure('Name', 'Computed allocation outputs from logged inputs');
    tiledlayout(7, 1, 'TileSpacing', 'compact');

    for channel = 1:7
        nexttile;
        hold on;
        legend_entries = {};

        if nargin >= 4 && ~isempty(output_log)
            if channel == 1
                plot(output_log.t_motor, output_log.motor0, 'k', 'LineWidth', 1.2);
            else
                plot(output_log.t_servo, output_log.servos(channel - 1, :), 'k', 'LineWidth', 1.2);
            end

            legend_entries{end+1} = "u_log"; %#ok<AGROW>
        end

        for case_idx = 1:numel(results)
            plot(t, results(case_idx).together.u(channel, :), '-.', 'LineWidth', 1.0);
            legend_entries{end+1} = results(case_idx).name + " together"; %#ok<AGROW>

            plot(t, results(case_idx).force_first.u(channel, :), '--', 'LineWidth', 1.0);
            legend_entries{end+1} = results(case_idx).name + " force-first"; %#ok<AGROW>
        end

        grid on;
        ylabel(channel_labels{channel});
        if channel == 1
            title('Computed u with logged actuator output overlay');
            legend(legend_entries, 'Interpreter', 'none', 'Location', 'eastoutside');
        end
    end
    xlabel('time [s]');
end

function plot_allocated_controls(t, c, results, axis_labels)
    figure('Name', 'Allocated control check');
    tiledlayout(6, 1, 'TileSpacing', 'compact');

    for axis_idx = 1:6
        nexttile;
        hold on;
        plot(t, c(axis_idx, :), 'k', 'LineWidth', 1.2);
        legend_entries = {"input"};

        for case_idx = 1:numel(results)
            plot(t, results(case_idx).together.allocated(axis_idx, :), 'LineWidth', 1.0);
            legend_entries{end+1} = results(case_idx).name + " together"; %#ok<AGROW>
            plot(t, results(case_idx).force_first.allocated(axis_idx, :), '--', 'LineWidth', 1.0);
            legend_entries{end+1} = results(case_idx).name + " force-first"; %#ok<AGROW>
        end

        grid on;
        ylabel(axis_labels{axis_idx});
        if axis_idx == 1
            title('allocated_control in PX4 scaled coordinates');
            legend(legend_entries, 'Interpreter', 'none', 'Location', 'eastoutside');
        end
    end
    xlabel('input time [s]');
end
