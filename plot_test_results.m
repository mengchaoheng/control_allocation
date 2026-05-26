function plot_test_results(results_file)
% Plot saved allocator comparison results without rerunning allocation.
%
% Usage:
%   test                  % computes and saves test_results.mat
%   plot_test_results     % redraws plots requested by reports_to_run
%   plot_test_results('test_results.mat')

if nargin < 1
    results_file = 'test_results.mat';
end

data = load(results_file);
if isfield(data, 'results')
    results = data.results;
else
    results = {data.result4, data.result6};
end

if isfield(data, 'reports_to_run')
    reports_to_run = data.reports_to_run;
else
    reports_to_run = {};
end

if isfield(data, 'use_restoring')
    plot_modes = get_plot_restoring_modes(data.use_restoring);
else
    plot_modes = {'restored'};
end

if has_report(reports_to_run, 'case_diff')
    if isfield(data, 'case_comparison_pairs') && isfield(data, 'allocation_methods_to_run')
        for mode_idx = 1:numel(plot_modes)
            plot_case_diff_u_curves(results, data.case_comparison_pairs, data.allocation_methods_to_run, data.t, data.len_command_px4, plot_modes{mode_idx});
        end
    end
end

if has_report(reports_to_run, 'method_diff')
    if isfield(data, 'allocation_methods_to_run')
        for mode_idx = 1:numel(plot_modes)
            plot_method_diff_u_curves(results, data.allocation_methods_to_run, data.t, data.len_command_px4, plot_modes{mode_idx});
        end
    end
end
end

function plot_case(result, t, len_command_px4)
    plot_algorithm(result.name, 'pca_dir vs cpp_dir', result, 'pca_dir', 'cpp_dir', t, len_command_px4);
    plot_algorithm(result.name, 'pca_dpscaled vs cpp_dpscaled', result, 'pca_dpscaled', 'cpp_dpscaled', t, len_command_px4);
    plot_algorithm(result.name, 'pca_prio vs cpp_prio', result, 'pca_prio', 'cpp_prio', t, len_command_px4);

    default_fields = {'name', 'cpp_tag', 'B', 'umin', 'umax', ...
        'pca_dir', 'pca_dpscaled', 'pca_prio', ...
        'pca_dir_raw', 'pca_dpscaled_raw', 'pca_prio_raw', ...
        'cpp_dir', 'cpp_dpscaled', 'cpp_prio', ...
        'cpp_dir_raw', 'cpp_dpscaled_raw', 'cpp_prio_raw'};
    fields = fieldnames(result);
    for idx = 1:numel(fields)
        field = fields{idx};
        if any(strcmp(field, default_fields))
            continue;
        end
        value = result.(field);
        if isnumeric(value) && ismatrix(value) && size(value, 1) == size(result.B, 2)
            plot_algorithm(result.name, field, result, field, '', t, len_command_px4);
        end
    end
end

function plot_algorithm(case_name, algorithm_name, result, matlab_field, cpp_field, t, len_command_px4)
    if ~isfield(result, matlab_field)
        return;
    end

    x_matlab = result.(matlab_field);
    n = min([len_command_px4, numel(t), size(x_matlab, 2)]);
    x_matlab = x_matlab(:, 1:n);
    has_cpp = ~isempty(cpp_field) && isfield(result, cpp_field);
    if has_cpp
        x_cpp = result.(cpp_field);
        n = min(n, size(x_cpp, 2));
        x_matlab = x_matlab(:, 1:n);
        x_cpp = x_cpp(:, 1:n);
    end

    figure('Name', [case_name ' ' algorithm_name]);
    tiledlayout(size(x_matlab, 1), 1, 'TileSpacing', 'compact');
    for idx = 1:size(x_matlab, 1)
        nexttile;
        if has_cpp
            plot(t(1:n), x_cpp(idx, :), 'r-', 'LineWidth', 1.1); hold on;
            plot(t(1:n), x_matlab(idx, :), 'b--', 'LineWidth', 1.0);
            legend('cpp', 'matlab', 'Location', 'best');
        else
            plot(t(1:n), x_matlab(idx, :), 'b-', 'LineWidth', 1.0);
            legend('matlab', 'Location', 'best');
        end
        ylabel(sprintf('u_%d', idx));
        grid on;
        if idx == 1
            title([case_name ' ' algorithm_name], 'Interpreter', 'none');
        end
        if idx == size(x_matlab, 1)
            xlabel('time [s]');
        end
    end
end

function plot_case_diff_u_curves(results, case_comparison_pairs, methods_to_run, t, len_command_px4, plot_mode)
    if isempty(case_comparison_pairs)
        return;
    end
    mode_label = get_plot_mode_label(plot_mode);

    for pair_idx = 1:size(case_comparison_pairs, 1)
        result_a = find_result_by_name(results, case_comparison_pairs{pair_idx, 1});
        result_b = find_result_by_name(results, case_comparison_pairs{pair_idx, 2});

        if isempty(result_a) || isempty(result_b)
            continue;
        end

        for method_idx = 1:numel(methods_to_run)
            method = methods_to_run{method_idx};
            u_a = get_case_method_u_for_plot(result_a, method, plot_mode);
            u_b = get_case_method_u_for_plot(result_b, method, plot_mode);

            if isempty(u_a) || isempty(u_b)
                continue;
            end

            n_outputs = max(size(u_a, 1), size(u_b, 1));
            u_a = pad_rows(u_a, n_outputs);
            u_b = pad_rows(u_b, n_outputs);
            n = min([len_command_px4, numel(t), size(u_a, 2), size(u_b, 2)]);

            figure('Name', sprintf('%s %s vs %s (%s)', method, result_a.name, result_b.name, mode_label));
            tiledlayout(n_outputs, 1, 'TileSpacing', 'compact');

            for output_idx = 1:n_outputs
                nexttile;
                styles = make_curve_styles({result_a.name, result_b.name});
                plot(t(1:n), u_a(output_idx, 1:n), ...
                     'Color', styles(1).color, ...
                     'LineStyle', styles(1).line_style, ...
                     'LineWidth', styles(1).line_width); hold on;
                plot(t(1:n), u_b(output_idx, 1:n), ...
                     'Color', styles(2).color, ...
                     'LineStyle', styles(2).line_style, ...
                     'LineWidth', styles(2).line_width);
                ylabel(sprintf('u_%d', output_idx));
                grid on;

                if output_idx == 1
                    title(sprintf('same allocator, different B: %s (%s)', method, mode_label), 'Interpreter', 'none');
                    legend(result_a.name, result_b.name, 'Location', 'best', 'Interpreter', 'none');
                end

                if output_idx == n_outputs
                    xlabel('time [s]');
                end
            end
        end
    end
end

function plot_method_diff_u_curves(results, methods_to_run, t, len_command_px4, plot_mode)
    mode_label = get_plot_mode_label(plot_mode);

    for case_idx = 1:numel(results)
        result = results{case_idx};
        method_curves = {};
        method_names = {};
        n_outputs = 0;

        for method_idx = 1:numel(methods_to_run)
            method = methods_to_run{method_idx};
            u = get_case_method_u_for_plot(result, method, plot_mode);

            if isempty(u)
                continue;
            end

            method_curves{end + 1} = u; %#ok<AGROW>
            method_names{end + 1} = method; %#ok<AGROW>
            n_outputs = max(n_outputs, size(u, 1));
        end

        if isempty(method_curves)
            continue;
        end

        figure('Name', sprintf('%s allocator u comparison (%s)', result.name, mode_label));
        tiledlayout(n_outputs, 1, 'TileSpacing', 'compact');
        line_styles = make_curve_styles(method_names);

        for output_idx = 1:n_outputs
            nexttile;
            hold on;

            for method_idx = 1:numel(method_curves)
                u = pad_rows(method_curves{method_idx}, n_outputs);
                n = min([len_command_px4, numel(t), size(u, 2)]);
                style = line_styles(method_idx);
                plot(t(1:n), u(output_idx, 1:n), ...
                     'Color', style.color, ...
                     'LineStyle', style.line_style, ...
                     'LineWidth', style.line_width);
            end

            ylabel(sprintf('u_%d', output_idx));
            grid on;

            if output_idx == 1
                title(sprintf('same B, different allocators: %s (%s)', result.name, mode_label), 'Interpreter', 'none');
                legend(method_names, 'Location', 'best', 'Interpreter', 'none');
            end

            if output_idx == n_outputs
                xlabel('time [s]');
            end
        end
    end
end

function styles = make_curve_styles(curve_names)
    % Stable per-curve styles for visual comparison.  Width stays in
    % [0.5, 1.5], while color and line style are also varied so two curves
    % in the same axes are easy to distinguish.
    line_styles = {'-', '--', ':', '-.'};
    color_order = [
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
        0.4660 0.6740 0.1880
        0.3010 0.7450 0.9330
        0.6350 0.0780 0.1840
        0.2500 0.2500 0.2500
    ];
    styles = struct('color', {}, 'line_style', {}, 'line_width', {});

    for curve_idx = 1:numel(curve_names)
        seed = sum(double(curve_names{curve_idx})) + 37 * curve_idx;
        width_random = mod(sin(seed * 12.9898) * 43758.5453, 1);

        if width_random < 0
            width_random = width_random + 1;
        end

        styles(curve_idx).color = color_order(mod(curve_idx - 1, size(color_order, 1)) + 1, :);
        styles(curve_idx).line_style = line_styles{mod(curve_idx - 1, numel(line_styles)) + 1};
        styles(curve_idx).line_width = 0.5 + width_random;
    end
end

function result = find_result_by_name(results, name)
    result = [];

    for idx = 1:numel(results)
        if strcmp(results{idx}.name, name)
            result = results{idx};
            return;
        end
    end
end

function u = get_case_method_u_for_plot(result, method, plot_mode)
    if nargin < 3
        plot_mode = 'restored';
    end

    use_raw = strcmpi(plot_mode, 'raw');

    if isfield(result, 'alloc') && isfield(result.alloc, method) && isfield(result.alloc.(method), 'u')
        if isfield(result.alloc.(method), 'ok') && ~any(result.alloc.(method).ok)
            u = [];
            return;
        end
        if use_raw && isfield(result.alloc.(method), 'u_raw')
            u = result.alloc.(method).u_raw;
        else
            u = result.alloc.(method).u;
        end
    elseif use_raw && isfield(result, [method '_raw'])
        u = result.([method '_raw']);
    elseif isfield(result, method)
        u = result.(method);
    else
        u = [];
        return;
    end

    if isfield(result, 'plot_num_outputs') && isfield(result, 'plot_active_idx') && numel(result.plot_active_idx) == size(u, 1)
        u_embedded = zeros(result.plot_num_outputs, size(u, 2));
        u_embedded(result.plot_active_idx, :) = u;
        u = u_embedded;
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
        label = 'no restoring';
    else
        label = 'restoring';
    end
end

function u = pad_rows(u, n_outputs)
    if size(u, 1) < n_outputs
        u(end + 1:n_outputs, :) = 0;
    end
end

function enabled = has_report(reports_to_run, report_name)
    if ischar(reports_to_run) || isstring(reports_to_run)
        reports_to_run = cellstr(reports_to_run);
    end

    enabled = any(strcmpi(reports_to_run, 'all')) || any(strcmpi(reports_to_run, report_name));
end
