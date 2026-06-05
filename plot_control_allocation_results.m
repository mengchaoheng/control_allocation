%% Plot control-allocation offline benchmark results
% 这个文件只负责画图，不重新解析日志，也不重新运行分配算法。
%
% 使用方式：
%   1. 先运行 main_control_allocation_benchmark.m，生成 results/*_allocation_compare_results.mat。
%   2. 再运行本文件：
%        plot_control_allocation_results
%
% 图的含义：
%   1. control_tracking:
%      同一段真实飞行指令 v_sp(t)，不同算法分配后得到的 B*u(t)。
%      如果 B*u 跟 v_sp 贴合，说明这组 B 和算法能满足该日志里的控制指令。
%      v_sp 固定顺序为 [Mx My Mz Fx Fy Fz]，由日志里的
%      vehicle_torque_setpoint + vehicle_thrust_setpoint 构造。
%
%   2. residual_v_minus_Bu:
%      residual = v_sp - B*u。
%      这是控制分配误差的主图，直接回答“算法分配出来以后还差多少”。
%
%   3. actuator_outputs:
%      每个执行器输出。PLOT_NORMALIZED=true 时画 u；
%      PLOT_NORMALIZED=false 时画 u * u_abs_max。
%      这个图用来观察不同算法在冗余执行器空间里怎么分配。
%      只有 u_px4 列数和当前 B 完全一致时，才叠加 PX4 日志输出。
%
%   4. summary:
%      residual 和 u 差异按 PLOT_NORMALIZED 的显示口径统计。
%      时间只看 allocator 平均耗时和 restoring 平均耗时，不画总运行时间；
%      分配误差画 mean/rms/max；离线 inv 参考画 rms/max；
%      fail/fallback 也直接画出来，方便定位算法异常点。
close all;
%% 用户可改配置
% RESULT_MAT 为空时，自动选择 RESULT_DIR 目录下最新的结果文件。
tool_dir = fileparts(mfilename('fullpath'));

if ~exist('RESULT_MAT', 'var')
    RESULT_MAT = "";
else
    RESULT_MAT = string(RESULT_MAT);
end

% METHODS_TO_PLOT 由 main 传进来：
%   METHODS_TO_PLOT = [METHODS_TO_RUN, CPP_ALLOCATOR_TO_PLOT]
% 空时画结果文件里已有的全部方法。
if ~exist('METHODS_TO_PLOT', 'var')
    METHODS_TO_PLOT = strings(0, 1);
end

% 执行器最多画 16 列，对应 PX4 ActuatorEffectiveness::NUM_ACTUATORS。
if ~exist('MAX_ACTUATORS_TO_PLOT', 'var')
    MAX_ACTUATORS_TO_PLOT = 16;
end

% 默认保存 png，同时打开图窗给你看。
% 如果你只想后台出图，把 SHOW_FIGURES 改成 false。
% 如果你还想保存 MATLAB 可交互修改的 .fig，把 SAVE_FIG 改成 true。
if ~exist('SAVE_PNG', 'var')
    SAVE_PNG = true;
end

if ~exist('SAVE_FIG', 'var')
    SAVE_FIG = false;
end

if ~exist('SHOW_FIGURES', 'var')
    SHOW_FIGURES = true;
end

% true 画归一化量；false 画能还原的非归一化量。
if ~exist('PLOT_NORMALIZED', 'var')
    PLOT_NORMALIZED = [];
end

%% 读取结果文件
% 这里不调用 main，也不调用 parser，保证画图不会把 benchmark 时间混进去。
if ~exist('RESULT_DIR', 'var') || strlength(string(RESULT_DIR)) == 0
    RESULT_DIR = fullfile(tool_dir, 'results');
end

if strlength(RESULT_MAT) == 0
    RESULT_MAT = select_latest_result_mat(RESULT_DIR);
end

if ~isfile(RESULT_MAT)
    error('找不到结果文件: %s。请先运行 main_control_allocation_benchmark。', RESULT_MAT);
end

S = load(RESULT_MAT);

results = S.results;
flightData = S.flightData;

if isfield(S, 'meta')
    meta = S.meta;
else
    meta = struct();
end

if isempty(PLOT_NORMALIZED)
    if isfield(meta, 'PLOT_NORMALIZED')
        PLOT_NORMALIZED = meta.PLOT_NORMALIZED;
    else
        PLOT_NORMALIZED = true;
    end
end

if isempty(METHODS_TO_PLOT) && isfield(meta, 'METHODS_TO_PLOT')
    METHODS_TO_PLOT = meta.METHODS_TO_PLOT;
end

[~, result_stem] = fileparts(char(RESULT_MAT));
figure_stem = regexprep(result_stem, '_allocation_compare_results$', '');

if ~exist('FIGURE_DIR', 'var') || strlength(string(FIGURE_DIR)) == 0
    FIGURE_DIR = fullfile(RESULT_DIR, [figure_stem '_figures']);
end

output_root = char(FIGURE_DIR);

if ~isfolder(output_root)
    mkdir(output_root);
end

fprintf('\nPlotting control-allocation result:\n  %s\n', RESULT_MAT);
fprintf('Figures will be saved to:\n  %s\n', output_root);
fprintf('Display units: %s\n', display_unit_text(PLOT_NORMALIZED));

if isfield(meta, 'benchmark_elapsed_s')
    fprintf('Benchmark wall time in result file: %.4f s\n', meta.benchmark_elapsed_s);
elseif isfield(S, 'benchmark_elapsed_s')
    fprintf('Benchmark wall time in result file: %.4f s\n', S.benchmark_elapsed_s);
end

if isfield(meta, 'PX4_REFERENCE_COMMIT')
    fprintf('PX4 reference commit: %s\n', meta.PX4_REFERENCE_COMMIT);
elseif isfield(S, 'PX4_REFERENCE_COMMIT')
    fprintf('PX4 reference commit: %s\n', S.PX4_REFERENCE_COMMIT);
end

%% 每个 B case 单独画一组图
% 一个 B case 对应一种 EffectivenessMatrix 输入。
% 未来如果同一机型有两个 allocator instance，就会自然形成两个 matrix/case，
% 或者一个 case 里包含不同维度的 u。这个画图文件只依赖结果结构，不假设机型。
saved_files = strings(0, 1);

for b_idx = 1:numel(results)
    case_result = results(b_idx);

    if ~isfield(case_result, 'alg') || isempty(case_result.alg)
        fprintf('Skip B case %d: no algorithm result.\n', b_idx);
        continue;
    end

    method_idx = select_method_indices(case_result.alg, METHODS_TO_PLOT);

    if isempty(method_idx)
        fprintf('Skip B case %s: no method selected.\n', case_result.name);
        continue;
    end

    t_rel = relative_time_s(case_result, flightData);
    active_axes = find_active_axes(case_result, method_idx, PLOT_NORMALIZED, meta);
    case_name = safe_file_name(string(case_result.name));
    case_output_dir = fullfile(output_root, char(case_name));

    if ~isfolder(case_output_dir)
        mkdir(case_output_dir);
    end

    fprintf('\nB case: %s\n', case_result.name);
    fprintf('  methods: %s\n', strjoin(string({case_result.alg(method_idx).name}), ', '));
    fprintf('  active axes: %s\n', mat2str(active_axes));

    fig = plot_control_tracking(case_result, method_idx, active_axes, t_rel, ...
        PLOT_NORMALIZED, meta, SHOW_FIGURES);
    saved_files = [saved_files; save_current_figure(fig, case_output_dir, ...
        case_name + "_control_tracking", SAVE_PNG, SAVE_FIG)]; %#ok<AGROW>

    fig = plot_residual_comparison(case_result, method_idx, active_axes, t_rel, ...
        PLOT_NORMALIZED, meta, SHOW_FIGURES);
    saved_files = [saved_files; save_current_figure(fig, case_output_dir, ...
        case_name + "_residual_v_minus_Bu", SAVE_PNG, SAVE_FIG)]; %#ok<AGROW>

    fig = plot_actuator_outputs(case_result, flightData, method_idx, t_rel, ...
        MAX_ACTUATORS_TO_PLOT, PLOT_NORMALIZED, SHOW_FIGURES);
    saved_files = [saved_files; save_current_figure(fig, case_output_dir, ...
        case_name + "_actuator_outputs", SAVE_PNG, SAVE_FIG)]; %#ok<AGROW>

    fig = plot_summary_metrics(case_result, method_idx, PLOT_NORMALIZED, meta, SHOW_FIGURES);
    saved_files = [saved_files; save_current_figure(fig, case_output_dir, ...
        case_name + "_summary", SAVE_PNG, SAVE_FIG)]; %#ok<AGROW>
end

comparison_figs = plot_same_algorithm_different_B_outputs(results, METHODS_TO_PLOT, ...
    MAX_ACTUATORS_TO_PLOT, PLOT_NORMALIZED, SHOW_FIGURES);

if ~isempty(comparison_figs)
    comparison_dir = fullfile(output_root, '_same_algorithm_different_B');

    if ~isfolder(comparison_dir)
        mkdir(comparison_dir);
    end

    for i = 1:numel(comparison_figs)
        saved_files = [saved_files; save_current_figure(comparison_figs(i).fig, ...
            comparison_dir, comparison_figs(i).stem, SAVE_PNG, SAVE_FIG)]; %#ok<AGROW>
    end
end

fprintf('\nSaved figures:\n');
for i = 1:numel(saved_files)
    fprintf('  %s\n', saved_files(i));
end

fprintf('\nDone.\n');

%% Local functions
function result_mat = select_latest_result_mat(result_dir)
% 自动选择最新结果文件，避免每次换日志都手动改路径。
files = dir(fullfile(result_dir, '*_allocation_compare_results.mat'));

if isempty(files)
    error('目录里没有 *_allocation_compare_results.mat: %s', result_dir);
end

[~, newest_idx] = max([files.datenum]);
result_mat = string(fullfile(files(newest_idx).folder, files(newest_idx).name));
end

function method_idx = select_method_indices(alg, methods_to_plot)
% 根据算法名筛选要画的算法。
names = string({alg.name});
methods_to_plot = string(methods_to_plot);
methods_to_plot = methods_to_plot(:).';

if isempty(methods_to_plot)
    method_idx = find(strlength(names) > 0);
else
    method_idx = find(ismember(names, methods_to_plot));
end
end

function t_rel = relative_time_s(case_result, flightData)
% 结果里优先用 case_result.t；如果旧结果没有这个字段，就退回 flightData.t。
if isfield(case_result, 't') && ~isempty(case_result.t)
    t = case_result.t(:);
elseif isfield(flightData, 't') && ~isempty(flightData.t)
    t = flightData.t(:);
elseif isfield(flightData, 'timestamp_us') && ~isempty(flightData.timestamp_us)
    t = double(flightData.timestamp_us(:)) * 1e-6;
else
    sample_count = size(case_result.alg(1).u, 1);
    t = (0:(sample_count - 1))';
end

t_rel = t - t(1);
end

function active_axes = find_active_axes(case_result, method_idx, plot_normalized, meta)
% 自动找有效控制轴。
% v = B*u 固定 6 轴：[Mx My Mz Fx Fy Fz]。
% 对 df4/SHC09 这类 surface allocation，通常只有前 3 个 torque 轴有量。
n_axes = 6;
axis_has_signal = false(1, n_axes);
v_sp = display_v_sp(case_result, plot_normalized, meta);

for i = method_idx(:)'
    alg = case_result.alg(i);
    v_achieved = display_v_achieved(case_result, alg, plot_normalized);
    residual = v_sp - v_achieved;
    data = abs(v_sp) + abs(v_achieved) + abs(residual);
    axis_has_signal = axis_has_signal | any(isfinite(data) & data > 1e-9, 1);
end

active_axes = find(axis_has_signal);

if isempty(active_axes)
    active_axes = 1:3;
end
end

function fig = plot_control_tracking(case_result, method_idx, active_axes, t_rel, ...
    plot_normalized, meta, show_figures)
% 画 v_sp 和 B*u 的时间曲线。
% 黑色虚线是同一条真实飞行指令；彩色实线是不同算法的实际分配效果。
fig = new_figure(show_figures, sprintf('%s control tracking', case_result.name));
layout = tiledlayout(numel(active_axes), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
unit_text = display_unit_text(plot_normalized);
title(layout, sprintf('%s: v_sp vs B*u (%s)', case_result.name, unit_text), 'Interpreter', 'none');
colors = lines(numel(method_idx));
labels = axis_labels();
v_sp = display_v_sp(case_result, plot_normalized, meta);

for tile_idx = 1:numel(active_axes)
    axis_id = active_axes(tile_idx);
    nexttile;
    hold on;
    grid on;

    plot(t_rel, v_sp(:, axis_id), 'k--', 'LineWidth', 1.2, 'DisplayName', 'v_sp');

    for k = 1:numel(method_idx)
        alg = case_result.alg(method_idx(k));
        v_achieved = display_v_achieved(case_result, alg, plot_normalized);
        plot(t_rel, v_achieved(:, axis_id), 'LineWidth', 1.0, ...
            'Color', colors(k, :), 'DisplayName', alg.name);
    end

    ylabel(sprintf('%s [%s]', labels(axis_id), unit_text), 'Interpreter', 'none');

    if tile_idx == numel(active_axes)
        xlabel('time since window start [s]');
    end

    if tile_idx == 1
        legend('Location', 'best', 'Interpreter', 'none');
    end
end
end

function fig = plot_residual_comparison(case_result, method_idx, active_axes, t_rel, ...
    plot_normalized, meta, show_figures)
% 画 residual = v_sp - B*u。
% 这个图是算法可行性和饱和情况最直接的观测。
fig = new_figure(show_figures, sprintf('%s residual', case_result.name));
layout = tiledlayout(numel(active_axes), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
unit_text = display_unit_text(plot_normalized);
title(layout, sprintf('%s: residual = v_sp - B*u (%s)', case_result.name, unit_text), 'Interpreter', 'none');
colors = lines(numel(method_idx));
labels = axis_labels();
v_sp = display_v_sp(case_result, plot_normalized, meta);

for tile_idx = 1:numel(active_axes)
    axis_id = active_axes(tile_idx);
    nexttile;
    hold on;
    grid on;

    for k = 1:numel(method_idx)
        alg = case_result.alg(method_idx(k));
        residual = v_sp - display_v_achieved(case_result, alg, plot_normalized);
        plot(t_rel, residual(:, axis_id), 'LineWidth', 1.0, ...
            'Color', colors(k, :), 'DisplayName', alg.name);
    end

    ylabel(sprintf('%s [%s]', labels(axis_id), unit_text), 'Interpreter', 'none');

    if tile_idx == numel(active_axes)
        xlabel('time since window start [s]');
    end

    if tile_idx == 1
        legend('Location', 'best', 'Interpreter', 'none');
    end
end
end

function fig = plot_actuator_outputs(case_result, flightData, method_idx, t_rel, ...
    max_actuators, plot_normalized, show_figures)
% 画执行器输出 u。
% u 是分配器归一化执行器空间里的 actuator_delta。
% 如果 PX4 日志输出的 actuator 列数和当前 B 完全一致，就用黑色虚线叠加为附加参考。
% 列数不一致通常表示日志机型/actuator 定义不同，此时只画离线算法结果。
first_alg = case_result.alg(method_idx(1));
u_dim = min([size(first_alg.u, 2), max_actuators]);
tile_rows = ceil(sqrt(u_dim));
tile_cols = ceil(u_dim / tile_rows);
colors = lines(numel(method_idx));
[u_ref, u_ref_name] = actuator_reference(case_result, flightData, u_dim);
u_scale = actuator_display_scale(case_result, u_dim, plot_normalized);
u_label = actuator_unit_text(plot_normalized);

fig = new_figure(show_figures, sprintf('%s actuator outputs', case_result.name));
layout = tiledlayout(tile_rows, tile_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('%s: actuator output %s', case_result.name, u_label), 'Interpreter', 'none');

for actuator = 1:u_dim
    nexttile;
    hold on;
    grid on;

    if ~isempty(u_ref)
        plot(t_rel, u_ref(:, actuator) * u_scale(actuator), 'k--', 'LineWidth', 1.0, ...
            'DisplayName', u_ref_name);
    end

    for k = 1:numel(method_idx)
        alg = case_result.alg(method_idx(k));
        plot(t_rel, alg.u(:, actuator) * u_scale(actuator), 'LineWidth', 0.9, ...
            'Color', colors(k, :), 'DisplayName', alg.name);
    end

    title(sprintf('u%d', actuator), 'Interpreter', 'none');

    if actuator > (tile_rows - 1) * tile_cols
        xlabel('time [s]');
    end

    ylabel(u_label);

    if actuator == 1
        legend('Location', 'best', 'Interpreter', 'none');
    end
end
end

function fig = plot_summary_metrics(case_result, method_idx, plot_normalized, meta, show_figures)
% 汇总图的数值口径和 print_control_allocation_comparison.m 保持一致。
% 这里不画 total elapsed time，因为它包含 MATLAB 绘图/函数调度/CSV 等非算法开销。
% 真正用于比较分配算法运行时间的是 avg_alloc(us/sample) 和 restore(us/sample)。
names = string({case_result.alg(method_idx).name});
avg_us = nan(size(method_idx));
restore_us = nan(size(method_idx));
residual_mean_abs = nan(size(method_idx));
residual_rms = nan(size(method_idx));
residual_max = nan(size(method_idx));
u_rmse_inv = nan(size(method_idx));
u_max_inv = nan(size(method_idx));
fail_count = nan(size(method_idx));
fallback_count = nan(size(method_idx));
sample_count = max(size(case_result.v_sp, 1), 1);
v_sp = display_v_sp(case_result, plot_normalized, meta);
u_scale = actuator_display_scale(case_result, size(case_result.B, 2), plot_normalized);

for k = 1:numel(method_idx)
    alg = case_result.alg(method_idx(k));
    avg_us(k) = alg.avg_us_per_sample;
    restore_us(k) = 1e6 * alg.restore_s / sample_count;
    residual = v_sp - display_v_achieved(case_result, alg, plot_normalized);
    residual_values = residual(isfinite(residual));

    if ~isempty(residual_values)
        residual_mean_abs(k) = mean(abs(residual_values));
    end

    [residual_rms(k), residual_max(k)] = finite_rmse_max(residual);
    if ~isempty(case_result.u_inv_ref)
        u_delta = (alg.u - case_result.u_inv_ref) .* u_scale;
        [u_rmse_inv(k), u_max_inv(k)] = finite_rmse_max(u_delta);
    else
        u_rmse_inv(k) = alg.rmse_vs_inv;
        u_max_inv(k) = alg.max_abs_vs_inv;
    end
    fail_count(k) = alg.fail_count;
    fallback_count(k) = alg.fallback_count;
end

fig = new_figure(show_figures, sprintf('%s summary', case_result.name));
layout = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
unit_text = display_unit_text(plot_normalized);
u_label = actuator_unit_text(plot_normalized);
title(layout, sprintf('%s: summary metrics (%s)', case_result.name, unit_text), 'Interpreter', 'none');

plot_metric_bar(names, avg_us, 'avg_alloc(us)', 'us / sample');
plot_metric_bar(names, restore_us, 'restore(us)', 'us / sample');
plot_metric_bar(names, residual_mean_abs, 'mean |v-Bu|', sprintf('mean abs [%s]', unit_text));
plot_metric_bar(names, residual_rms, 'rms(v-Bu)', sprintf('rms [%s]', unit_text));
plot_metric_bar(names, residual_max, 'max |v-Bu|', sprintf('max abs [%s]', unit_text));
plot_metric_bar(names, u_rmse_inv, 'rms(u-u_inv)', sprintf('rms [%s]', u_label));
plot_metric_bar(names, u_max_inv, 'max |u-u_inv|', sprintf('max abs [%s]', u_label));
plot_metric_bar(names, fail_count, 'fail', 'samples');
plot_metric_bar(names, fallback_count, 'fallback', 'samples');
end

function plot_metric_bar(names, values, title_text, ylabel_text)
% 一个小工具，保证 summary 的每个子图样式和标签一致。
nexttile;
bar(values);
grid on;
title(title_text, 'Interpreter', 'none');
ylabel(ylabel_text, 'Interpreter', 'none');

finite_values = values(isfinite(values));

if isempty(finite_values) || max(finite_values) <= 0
    ylim([0 1]);
else
    ylim([0 1.08 * max(finite_values)]);
end

set(gca, 'XTick', 1:numel(names), ...
    'XTickLabel', names, ...
    'XTickLabelRotation', 30, ...
    'TickLabelInterpreter', 'none');
end

function figures = plot_same_algorithm_different_B_outputs(results, methods_to_plot, ...
    max_actuators, plot_normalized, show_figures)
% 同一个算法，不同 B case 的 u 输出叠图。用于比较“换 B 以后同一算法怎么分配”。
figures = struct('fig', {}, 'stem', {});

if numel(results) < 2
    return;
end

common_names = string({results(1).alg.name});

for case_idx = 2:numel(results)
    common_names = intersect(common_names, string({results(case_idx).alg.name}), 'stable');
end

methods_to_plot = string(methods_to_plot);
methods_to_plot = methods_to_plot(:).';
if ~isempty(methods_to_plot)
    common_names = common_names(ismember(common_names, methods_to_plot));
end

for method_idx = 1:numel(common_names)
    method_name = common_names(method_idx);
    alg_idx = zeros(1, numel(results));
    u_dim = 0;

    for case_idx = 1:numel(results)
        names = string({results(case_idx).alg.name});
        alg_idx(case_idx) = find(names == method_name, 1);
        u_dim = max(u_dim, size(results(case_idx).alg(alg_idx(case_idx)).u, 2));
    end

    u_dim = min(u_dim, max_actuators);
    tile_rows = ceil(sqrt(u_dim));
    tile_cols = ceil(u_dim / tile_rows);
    colors = lines(numel(results));
    fig = new_figure(show_figures, sprintf('%s: same algorithm, different B', method_name));
    layout = tiledlayout(tile_rows, tile_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    u_label = actuator_unit_text(plot_normalized);
    title(layout, sprintf('same algorithm, different B: %s (%s)', method_name, u_label), 'Interpreter', 'none');

    for actuator = 1:u_dim
        nexttile;
        hold on;
        grid on;

        for case_idx = 1:numel(results)
            alg = results(case_idx).alg(alg_idx(case_idx));

            if size(alg.u, 2) < actuator
                continue;
            end

            t_rel = relative_time_s(results(case_idx), struct());
            u_scale = actuator_display_scale(results(case_idx), u_dim, plot_normalized);
            plot(t_rel, alg.u(:, actuator) * u_scale(actuator), 'LineWidth', 0.9, ...
                'Color', colors(case_idx, :), 'DisplayName', results(case_idx).name);
        end

        title(sprintf('u%d', actuator), 'Interpreter', 'none');
        ylabel(u_label);

        if actuator > (tile_rows - 1) * tile_cols
            xlabel('time [s]');
        end

        if actuator == 1
            legend('Location', 'best', 'Interpreter', 'none');
        end
    end

    figures(end+1).fig = fig; %#ok<AGROW>
    figures(end).stem = safe_file_name("same_algorithm_" + method_name + "_u_by_B");
end
end

function [u_ref, u_ref_name] = actuator_reference(case_result, flightData, u_dim)
% 找 PX4 在线执行器输出。parser 已经插值到 command 时间轴。
% 必须要求列数和当前 B 完全一致；否则 actuator 顺序/定义不一定对应。
u_ref = [];
u_ref_name = "";

if ~isempty(flightData.u_px4) ...
        && size(flightData.u_px4, 2) == size(case_result.B, 2) ...
        && size(flightData.u_px4, 2) >= u_dim
    u_ref = flightData.u_px4(:, 1:u_dim);
    u_ref_name = "u_px4";
end
end

function v_sp = display_v_sp(case_result, plot_normalized, meta)
v_sp = fixed_width(case_result.v_sp, 6);

if ~plot_normalized && isfield(meta, 'NORMALIZE_V_SP') && meta.NORMALIZE_V_SP
    D = case_D(case_result);
    v_sp = (D \ v_sp')';
end
end

function v_achieved = display_v_achieved(case_result, alg, plot_normalized)
v_achieved = fixed_width(alg.v_achieved, 6);

if ~plot_normalized
    D = case_D(case_result);
    v_achieved = (D \ v_achieved')';
end
end

function D = case_D(case_result)
if isfield(case_result, 'D') && isequal(size(case_result.D), [6 6])
    D = case_result.D;
else
    D = eye(6);
end
end

function u_scale = actuator_display_scale(case_result, u_dim, plot_normalized)
u_scale = ones(1, u_dim);

if plot_normalized || ~isfield(case_result, 'u_abs_max') || isempty(case_result.u_abs_max)
    return;
end

u_abs_max = double(case_result.u_abs_max(:))';
n = min(u_dim, numel(u_abs_max));
u_scale(1:n) = u_abs_max(1:n);
end

function unit_text = display_unit_text(plot_normalized)
if plot_normalized
    unit_text = 'normalized';
else
    unit_text = 'non-normalized';
end
end

function unit_text = actuator_unit_text(plot_normalized)
if plot_normalized
    unit_text = 'u';
else
    unit_text = 'u * u_abs_max';
end
end

function [rmse_value, max_value] = finite_rmse_max(X)
values = X(isfinite(X));

if isempty(values)
    rmse_value = nan;
    max_value = nan;
else
    rmse_value = sqrt(mean(values(:).^2));
    max_value = max(abs(values(:)));
end
end

function X = fixed_width(X, n_cols)
% 把矩阵补齐/裁剪到固定列数，方便 6 轴图统一处理。
if isempty(X)
    X = nan(0, n_cols);
    return;
end

if size(X, 2) < n_cols
    X = [X nan(size(X, 1), n_cols - size(X, 2))];
elseif size(X, 2) > n_cols
    X = X(:, 1:n_cols);
end
end

function labels = axis_labels()
% v = B*u 的 6 维顺序。
labels = ["Mx"; "My"; "Mz"; "Fx"; "Fy"; "Fz"];
end

function fig = new_figure(show_figures, fig_name)
% 默认 headless 保存图；需要交互查看时把 SHOW_FIGURES 改成 true。
if show_figures
    visible = 'on';
else
    visible = 'off';
end

fig = figure('Visible', visible, 'Name', fig_name, 'Color', 'w');
fig.Position(3:4) = [1400 900];
end

function saved_files = save_current_figure(fig, output_dir, file_stem, save_png, save_fig)
% 同一张图保存 png/fig 两种格式。
saved_files = strings(0, 1);

if save_png
    png_path = fullfile(output_dir, char(file_stem + ".png"));
    exportgraphics(fig, png_path, 'Resolution', 160);
    saved_files(end+1, 1) = string(png_path); %#ok<AGROW>
end

if save_fig
    fig_path = fullfile(output_dir, char(file_stem + ".fig"));
    savefig(fig, fig_path);
    saved_files(end+1, 1) = string(fig_path); %#ok<AGROW>
end

% close(fig);
end

function name = safe_file_name(name)
% 把 B case 名字变成可以作为文件名的字符串。
name = regexprep(name, '[^a-zA-Z0-9_\-]+', '_');
name = regexprep(name, '_+', '_');
name = strip(name, '_');

if strlength(name) == 0
    name = "case";
end
end
