function plot_test_results(results_file)
% Plot saved allocator comparison results without rerunning allocation.
%
% Usage:
%   test                  % computes and saves test_results.mat
%   plot_test_results     % redraws from the saved data
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

for case_idx = 1:numel(results)
    result = results{case_idx};
    plot_case(result, data.t, data.len_command_px4);
end
end

function plot_case(result, t, len_command_px4)
    plot_algorithm(result.name, 'DP_LPCA', result, 'matlab_dp', 'cpp_dp', t, len_command_px4);
    plot_algorithm(result.name, 'DPscaled_LPCA', result, 'matlab_dpscaled', 'cpp_dpscaled', t, len_command_px4);
    plot_algorithm(result.name, 'DP_LPCA_prio', result, 'matlab_prio', 'cpp_prio', t, len_command_px4);

    default_fields = {'name', 'cpp_tag', 'B', 'umin', 'umax', ...
        'matlab_dp', 'matlab_dpscaled', 'matlab_prio', ...
        'matlab_dp_raw', 'matlab_dpscaled_raw', 'matlab_prio_raw', ...
        'cpp_dp', 'cpp_dpscaled', 'cpp_prio'};
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
