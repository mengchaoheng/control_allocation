function plot_test_results(results_file)
%PLOT_TEST_RESULTS Plot allocator outputs and allocation errors saved by test.m.
%
% 这个文件放在 control_allocation 根目录，避免 test.m 继续调用
% _pending_delete_legacy_mfiles/plot_test_results.m 里的旧图。
%
% 图的口径：
%   1) u comparison:
%      同一个 B case 下，不同分配算法的执行器输出 u(t)。
%   2) tracking:
%      同一张图画 y(t) 和 B*u(t)，直接看跟踪情况。
%   3) allocation error:
%      err(t) = B*u(t) - y(t)，以及 ||err(t)||_2。
%      这是判断算法是否真正达到控制输入的核心图。

if nargin < 1
    results_file = 'test_results.mat';
end

close all;

data = load(results_file);
results = read_results_cell(data);

if ~isfield(data, 'command_px4') || ~isfield(data, 't')
    error('plot_test_results:MissingData', ...
          '%s must contain command_px4 and t saved by test.m.', results_file);
end

command = data.command_px4;
t = data.t;
methods_to_run = read_field_or_default(data, 'allocation_methods_to_run', {});
reports_to_run = read_field_or_default(data, 'reports_to_run', {'method_diff'});
plot_modes = get_plot_restoring_modes(read_field_or_default(data, 'use_restoring', 'restored'));

[results_dir, results_stem] = fileparts(results_file);
if isempty(results_dir)
    results_dir = pwd;
end
fig_dir = fullfile(results_dir, [results_stem '_figures']);
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

fprintf('\nPlotting test results from %s\n', results_file);
fprintf('  figure dir: %s\n', fig_dir);

% test.m 里 method_diff 是最常用报告；这里即使 reports_to_run 没写 method_diff，
% 只要有 methods_to_run，也画误差，因为分配误差是本测试最关键的图。
if isempty(methods_to_run)
    methods_to_run = infer_methods_from_results(results);
end

for mode_idx = 1:numel(plot_modes)
    plot_mode = plot_modes{mode_idx};

    for case_idx = 1:numel(results)
        result = results{case_idx};
        if ~isfield(result, 'alloc') || isempty(methods_to_run)
            continue;
        end

        plot_method_u_comparison(result, methods_to_run, t, command, plot_mode, fig_dir);
        plot_method_tracking(result, methods_to_run, t, command, plot_mode, fig_dir);
        plot_method_allocation_error(result, methods_to_run, t, command, plot_mode, fig_dir);
    end
end

if has_report(reports_to_run, 'case_diff') && isfield(data, 'case_comparison_pairs')
    for mode_idx = 1:numel(plot_modes)
        plot_case_error_comparison(results, data.case_comparison_pairs, methods_to_run, t, command, plot_modes{mode_idx}, fig_dir);
    end
end
end

function results = read_results_cell(data)
% test.m 新版保存 results；旧版可能只保存 result4/result6。
if isfield(data, 'results')
    results = data.results;
elseif isfield(data, 'result4') || isfield(data, 'result6')
    results = {};
    if isfield(data, 'result4') && ~isempty(data.result4)
        results{end + 1} = data.result4; %#ok<AGROW>
    end
    if isfield(data, 'result6') && ~isempty(data.result6)
        results{end + 1} = data.result6; %#ok<AGROW>
    end
else
    error('plot_test_results:MissingResults', 'No results found in MAT file.');
end
end

function value = read_field_or_default(data, field_name, default_value)
if isfield(data, field_name)
    value = data.(field_name);
else
    value = default_value;
end
end

function methods = infer_methods_from_results(results)
% 没有保存 allocation_methods_to_run 时，从 result.alloc 字段推断可画算法。
methods = {};
for case_idx = 1:numel(results)
    if isfield(results{case_idx}, 'alloc')
        names = fieldnames(results{case_idx}.alloc);
        methods = union(methods, names, 'stable');
    end
end
end

function plot_method_u_comparison(result, methods_to_run, t, command, plot_mode, fig_dir)
% 画执行器输出 u。它回答“两个算法输出差多少”，但不直接说明谁更好。
[curves, names] = collect_method_u_curves(result, methods_to_run, plot_mode);
if isempty(curves)
    return;
end

n_outputs = max(cellfun(@(u) size(u, 1), curves));
n = common_sample_count(t, command, curves);
styles = make_curve_styles(names);
mode_label = get_plot_mode_label(plot_mode);

fig = new_figure(sprintf('%s u comparison (%s)', result.name, mode_label));
tiledlayout(n_outputs, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for output_idx = 1:n_outputs
    nexttile;
    hold on;

    for method_idx = 1:numel(curves)
        u = pad_rows(curves{method_idx}, n_outputs);
        plot(t(1:n), u(output_idx, 1:n), ...
             'Color', styles(method_idx).color, ...
             'LineStyle', styles(method_idx).line_style, ...
             'LineWidth', styles(method_idx).line_width);
    end

    ylabel(sprintf('u_%d', output_idx));
    grid on;

    if output_idx == 1
        title(sprintf('%s: allocator output u (%s)', result.name, mode_label), 'Interpreter', 'none');
        legend(names, 'Location', 'best', 'Interpreter', 'none');
    end

    if output_idx == n_outputs
        xlabel('time [s]');
    end
end

save_or_show_figure(fig, fig_dir, sprintf('%s_u_%s', result.name, mode_label));
end

function plot_method_allocation_error(result, methods_to_run, t, command, plot_mode, fig_dir)
% 画 err = B*u - y。这个图才是分配算法是否满足输入的直接证据。
[curves, names] = collect_method_u_curves(result, methods_to_run, plot_mode);
if isempty(curves)
    return;
end

n = common_sample_count(t, command, curves);
B = result.B;
y = command(:, 1:n);
styles = make_curve_styles(names);
mode_label = get_plot_mode_label(plot_mode);
n_axes = size(y, 1);

fig = new_figure(sprintf('%s allocation error (%s)', result.name, mode_label));
tiledlayout(n_axes + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for axis_idx = 1:n_axes
    nexttile;
    hold on;

    for method_idx = 1:numel(curves)
        u = curves{method_idx}(:, 1:n);
        err = B * u - y;
        plot(t(1:n), err(axis_idx, :), ...
             'Color', styles(method_idx).color, ...
             'LineStyle', styles(method_idx).line_style, ...
             'LineWidth', styles(method_idx).line_width);
    end

    ylabel(sprintf('e_%d', axis_idx));
    grid on;

    if axis_idx == 1
        title(sprintf('%s: allocation error B*u - y (%s)', result.name, mode_label), 'Interpreter', 'none');
        legend(names, 'Location', 'best', 'Interpreter', 'none');
    end
end

nexttile;
hold on;
for method_idx = 1:numel(curves)
    u = curves{method_idx}(:, 1:n);
    err_norm = vecnorm(B * u - y, 2, 1);
    plot(t(1:n), err_norm, ...
         'Color', styles(method_idx).color, ...
         'LineStyle', styles(method_idx).line_style, ...
         'LineWidth', max(styles(method_idx).line_width, 1.1));
end
ylabel('||e||_2');
xlabel('time [s]');
grid on;

save_or_show_figure(fig, fig_dir, sprintf('%s_allocation_error_%s', result.name, mode_label));
print_allocation_error_summary(result, methods_to_run, t, command, plot_mode);
end

function plot_method_tracking(result, methods_to_run, t, command, plot_mode, fig_dir)
% 画 y command 和 B*u achieved。它回答“算法有没有跟上输入”。
[curves, names] = collect_method_u_curves(result, methods_to_run, plot_mode);
if isempty(curves)
    return;
end

n = common_sample_count(t, command, curves);
B = result.B;
y = command(:, 1:n);
styles = make_curve_styles(names);
mode_label = get_plot_mode_label(plot_mode);
n_axes = size(y, 1);

fig = new_figure(sprintf('%s tracking (%s)', result.name, mode_label));
tiledlayout(n_axes, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for axis_idx = 1:n_axes
    nexttile;
    hold on;

    % command 用黑线，所有算法的 achieved control 用不同颜色。
    plot(t(1:n), y(axis_idx, :), 'k-', 'LineWidth', 1.35);

    for method_idx = 1:numel(curves)
        u = curves{method_idx}(:, 1:n);
        achieved = B * u;
        plot(t(1:n), achieved(axis_idx, :), ...
             'Color', styles(method_idx).color, ...
             'LineStyle', styles(method_idx).line_style, ...
             'LineWidth', styles(method_idx).line_width);
    end

    ylabel(sprintf('y_%d', axis_idx));
    grid on;

    if axis_idx == 1
        title(sprintf('%s: command tracking, y vs B*u (%s)', result.name, mode_label), 'Interpreter', 'none');
        legend(['command', names], 'Location', 'best', 'Interpreter', 'none');
    end

    if axis_idx == n_axes
        xlabel('time [s]');
    end
end

save_or_show_figure(fig, fig_dir, sprintf('%s_tracking_%s', result.name, mode_label));
end

function plot_case_error_comparison(results, case_pairs, methods_to_run, t, command, plot_mode, fig_dir)
% 可选：同一个算法在不同 B case 下的误差范数对比。
mode_label = get_plot_mode_label(plot_mode);

for pair_idx = 1:size(case_pairs, 1)
    result_a = find_result_by_name(results, case_pairs{pair_idx, 1});
    result_b = find_result_by_name(results, case_pairs{pair_idx, 2});
    if isempty(result_a) || isempty(result_b)
        continue;
    end

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};
        u_a = get_method_u(result_a, method, plot_mode);
        u_b = get_method_u(result_b, method, plot_mode);
        if isempty(u_a) || isempty(u_b)
            continue;
        end

        n = min([numel(t), size(command, 2), size(u_a, 2), size(u_b, 2)]);
        err_a = vecnorm(result_a.B * u_a(:, 1:n) - command(:, 1:n), 2, 1);
        err_b = vecnorm(result_b.B * u_b(:, 1:n) - command(:, 1:n), 2, 1);

        fig = new_figure(sprintf('%s case error comparison (%s)', method, mode_label));
        plot(t(1:n), err_a, 'LineWidth', 1.2); hold on;
        plot(t(1:n), err_b, 'LineWidth', 1.2);
        grid on;
        xlabel('time [s]');
        ylabel('||B*u-y||_2');
        title(sprintf('%s: allocation error norm, different B (%s)', method, mode_label), 'Interpreter', 'none');
        legend({result_a.name, result_b.name}, 'Location', 'best', 'Interpreter', 'none');

        save_or_show_figure(fig, fig_dir, sprintf('%s_%s_vs_%s_error_%s', method, result_a.name, result_b.name, mode_label));
    end
end
end

function print_allocation_error_summary(result, methods_to_run, t, command, plot_mode)
% 图之外再打印一张短表，避免只靠肉眼读曲线。
mode_label = get_plot_mode_label(plot_mode);
fprintf('\nAllocation error summary: %s (%s)\n', result.name, mode_label);
fprintf('  %-14s %12s %12s %12s %12s %12s\n', ...
        'method', 'max|e|', 'rms||e||', 'max||e||', 'worst_t', 'max|u|');

for method_idx = 1:numel(methods_to_run)
    method = methods_to_run{method_idx};
    u = get_method_u(result, method, plot_mode);
    if isempty(u)
        continue;
    end

    n = min([numel(t), size(command, 2), size(u, 2)]);
    err = result.B * u(:, 1:n) - command(:, 1:n);
    err_norm = vecnorm(err, 2, 1);
    [max_norm, worst_idx] = max(err_norm);
    fprintf('  %-14s %12.6g %12.6g %12.6g %12.6g %12.6g\n', ...
            method, max(abs(err), [], 'all'), sqrt(mean(err_norm.^2)), ...
            max_norm, t(worst_idx), max(abs(u(:, 1:n)), [], 'all'));
end
end

function [curves, names] = collect_method_u_curves(result, methods_to_run, plot_mode)
curves = {};
names = {};

for method_idx = 1:numel(methods_to_run)
    method = methods_to_run{method_idx};
    u = get_method_u(result, method, plot_mode);
    if isempty(u)
        continue;
    end

    curves{end + 1} = u; %#ok<AGROW>
    names{end + 1} = char(method); %#ok<AGROW>
end
end

function u = get_method_u(result, method, plot_mode)
% raw 模式用 u_raw；restored 模式用 u。没有 ok 样本的算法不画。
method = char(method);
use_raw = strcmpi(plot_mode, 'raw');

if ~isfield(result, 'alloc') || ~isfield(result.alloc, method)
    u = [];
    return;
end

data = result.alloc.(method);
if isfield(data, 'ok') && ~any(data.ok)
    u = [];
    return;
end

if use_raw && isfield(data, 'u_raw')
    u = data.u_raw;
else
    u = data.u;
end
end

function n = common_sample_count(t, command, curves)
n = min([numel(t), size(command, 2), cellfun(@(u) size(u, 2), curves)]);
end

function fig = new_figure(name)
% batch/headless 运行时也保存 PNG；桌面 MATLAB 运行时同时弹图。
if usejava('desktop')
    visible = 'on';
else
    visible = 'off';
end
fig = figure('Name', name, 'Visible', visible);
end

function save_or_show_figure(fig, fig_dir, fig_name)
png_path = fullfile(fig_dir, [sanitize_filename(fig_name) '.png']);
exportgraphics(fig, png_path, 'Resolution', 180);
fprintf('  saved figure: %s\n', png_path);

if ~usejava('desktop')
    close(fig);
end
end

function safe = sanitize_filename(name)
safe = regexprep(char(name), '[^\w\-.]+', '_');
safe = regexprep(safe, '_+', '_');
safe = regexprep(safe, '^_|_$', '');
end

function styles = make_curve_styles(curve_names)
% 固定颜色/线型，保证每次重画同一算法看起来一致。
line_styles = {'-', '--', ':', '-.'};
color_order = [
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.4660 0.6740 0.1880
    0.4940 0.1840 0.5560
    0.9290 0.6940 0.1250
    0.3010 0.7450 0.9330
    0.6350 0.0780 0.1840
    0.2500 0.2500 0.2500
];
styles = struct('color', {}, 'line_style', {}, 'line_width', {});

for curve_idx = 1:numel(curve_names)
    styles(curve_idx).color = color_order(mod(curve_idx - 1, size(color_order, 1)) + 1, :);
    styles(curve_idx).line_style = line_styles{mod(curve_idx - 1, numel(line_styles)) + 1};
    styles(curve_idx).line_width = 1.05;
end
end

function result = find_result_by_name(results, name)
result = [];
for idx = 1:numel(results)
    if isfield(results{idx}, 'name') && strcmp(results{idx}.name, name)
        result = results{idx};
        return;
    end
end
end

function u = pad_rows(u, n_outputs)
if size(u, 1) < n_outputs
    u(end + 1:n_outputs, :) = 0;
end
end

function plot_modes = get_plot_restoring_modes(use_restoring)
if ischar(use_restoring) || isstring(use_restoring)
    mode = lower(char(use_restoring));
    if strcmp(mode, 'both')
        plot_modes = {'restored', 'raw'};
    elseif any(strcmp(mode, {'false', 'off', 'no', 'raw'}))
        plot_modes = {'raw'};
    else
        plot_modes = {'restored'};
    end
elseif use_restoring
    plot_modes = {'restored'};
else
    plot_modes = {'raw'};
end
end

function label = get_plot_mode_label(plot_mode)
if strcmpi(plot_mode, 'raw')
    label = 'raw';
else
    label = 'restored';
end
end

function enabled = has_report(reports_to_run, report_name)
if ischar(reports_to_run) || isstring(reports_to_run)
    reports_to_run = cellstr(reports_to_run);
end

enabled = any(strcmpi(reports_to_run, 'all')) || any(strcmpi(reports_to_run, report_name));
end
