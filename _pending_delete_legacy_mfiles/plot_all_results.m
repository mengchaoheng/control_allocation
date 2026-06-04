function plot_all_results(t, u_reference, results, channel_labels, B_eval, v_eval, axis_labels, reference_name, u_reference_label, cpp_results)
% 打开控制分配测试的所有关键图，不保存 PNG。
%   1) 每个执行器通道 u_reference 和各算法 u 的输出曲线。
%   2) 每个执行器通道 u_algorithm-u_reference 的误差曲线。
%   3) 每个控制轴 B*u-v 的 tracking error。
%   4) 每个执行器通道相对 reference_name 的输出差异。
%   5) C++ replay 和 MATLAB 同名算法的输出叠图与差异。
    if nargin < 9 || isempty(u_reference_label)
        u_reference_label = 'u_reference';
    end
    if nargin < 10 || isempty(cpp_results)
        cpp_results = struct();
    end

    plot_actuator_outputs(t, u_reference, results, channel_labels, u_reference_label);
    plot_actuator_reference_line_errors(t, u_reference, results, channel_labels, u_reference_label);
    plot_tracking_errors(t, results, B_eval, v_eval, axis_labels);
    plot_actuator_reference_errors(t, results, channel_labels, reference_name);
    plot_cpp_replay_actuator_overlays(t, results, cpp_results, channel_labels);
    plot_cpp_replay_actuator_differences(t, results, cpp_results, channel_labels);
    plot_cpp_replay_tracking_differences(t, results, cpp_results, B_eval, axis_labels);
end

function plot_actuator_outputs(t, u_reference, results, channel_labels, u_reference_label)
    names = fieldnames(results);
    figure('Name', 'SHC09 allocation outputs', 'Color', 'w');

    for ch = 1:numel(channel_labels)
        subplot(numel(channel_labels), 1, ch);
        plot(t, u_reference(ch, :), 'k-', 'LineWidth', 0.8, 'DisplayName', u_reference_label);
        hold on;

        for i = 1:numel(names)
            plot(t, results.(names{i}).u(ch, :), ...
                'LineWidth', 0.7, 'DisplayName', names{i});
        end

        grid on;
        ylabel(sprintf('u%d %s', ch, channel_labels{ch}), 'Interpreter', 'none');

        if ch == 1
            title('SHC09 控制分配输出对比', 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 3, 'Interpreter', 'none');
        end
        if ch == numel(channel_labels)
            xlabel('Time (s)');
        end
    end
end

function plot_actuator_reference_line_errors(t, u_reference, results, channel_labels, u_reference_label)
    names = fieldnames(results);
    figure('Name', ['SHC09 actuator error vs ' u_reference_label], 'Color', 'w');

    for ch = 1:numel(channel_labels)
        subplot(numel(channel_labels), 1, ch);
        hold on;

        for i = 1:numel(names)
            err = results.(names{i}).u(ch, :) - u_reference(ch, :);
            plot(t, err, 'LineWidth', 0.7, 'DisplayName', names{i});
        end

        yline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        grid on;
        ylabel(sprintf('du%d %s', ch, channel_labels{ch}), 'Interpreter', 'none');

        if ch == 1
            title(['执行器输出误差: u_algorithm - ' u_reference_label], 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 3, 'Interpreter', 'none');
        end
        if ch == numel(channel_labels)
            xlabel('Time (s)');
        end
    end
end

function plot_tracking_errors(t, results, B_eval, v_eval, axis_labels)
    names = fieldnames(results);
    figure('Name', 'SHC09 tracking error B*u-v', 'Color', 'w');

    for ax = 1:numel(axis_labels)
        subplot(numel(axis_labels), 1, ax);
        hold on;

        for i = 1:numel(names)
            u = results.(names{i}).u;
            err = B_eval(ax, :) * u - v_eval(ax, :);
            plot(t, err, 'LineWidth', 0.7, 'DisplayName', names{i});
        end

        yline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        grid on;
        ylabel(axis_labels{ax}, 'Interpreter', 'none');

        if ax == 1
            title('控制量跟踪误差: B_{par}*u - v_{par}', 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 3, 'Interpreter', 'none');
        end
        if ax == numel(axis_labels)
            xlabel('Time (s)');
        end
    end
end

function plot_actuator_reference_errors(t, results, channel_labels, reference_name)
    if isempty(reference_name) || ~isfield(results, reference_name)
        return;
    end

    names = fieldnames(results);
    u_ref = results.(reference_name).u;
    figure('Name', ['SHC09 actuator difference vs ' reference_name], 'Color', 'w');

    for ch = 1:numel(channel_labels)
        subplot(numel(channel_labels), 1, ch);
        hold on;

        for i = 1:numel(names)
            name = names{i};
            err = results.(name).u(ch, :) - u_ref(ch, :);
            plot(t, err, 'LineWidth', 0.7, 'DisplayName', name);
        end

        yline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        grid on;
        ylabel(sprintf('du%d %s', ch, channel_labels{ch}), 'Interpreter', 'none');

        if ch == 1
            title(['执行器输出差异: u - ' reference_name], 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 3, 'Interpreter', 'none');
        end
        if ch == numel(channel_labels)
            xlabel('Time (s)');
        end
    end
end

function plot_cpp_replay_actuator_overlays(t, matlab_results, cpp_results, channel_labels)
    names = common_result_names(matlab_results, cpp_results);
    if isempty(names)
        report_cpp_plot_skipped(matlab_results, cpp_results);
        return;
    end

    colors = lines(numel(names));
    figure('Name', 'SHC09 C++ replay actuator overlay', 'Color', 'w');

    for ch = 1:numel(channel_labels)
        subplot(numel(channel_labels), 1, ch);
        hold on;

        for i = 1:numel(names)
            name = names{i};
            n = common_sample_count(t, matlab_results.(name), cpp_results.(name));
            good = common_success_mask(matlab_results.(name), cpp_results.(name), n);
            matlab_u = matlab_results.(name).u(ch, 1:n);
            cpp_u = cpp_results.(name).u(ch, 1:n);
            matlab_u(~good) = NaN;
            cpp_u(~good) = NaN;

            plot(t(1:n), matlab_u, ...
                '-', 'Color', colors(i, :), 'LineWidth', 0.7, ...
                'DisplayName', [name ' matlab']);
            plot(t(1:n), cpp_u, ...
                '--', 'Color', colors(i, :), 'LineWidth', 0.7, ...
                'DisplayName', [name ' cpp']);

            if any(~good)
                plot(t(~good), cpp_results.(name).u(ch, ~good), ...
                    'x', 'Color', colors(i, :), 'LineWidth', 0.7, ...
                    'HandleVisibility', 'off');
            end
        end

        grid on;
        ylabel(sprintf('u%d %s', ch, channel_labels{ch}), 'Interpreter', 'none');

        if ch == 1
            title('C++ replay 对比: MATLAB 实线 / C++ 虚线', 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 2, 'Interpreter', 'none');
        end
        if ch == numel(channel_labels)
            xlabel('Time (s)');
        end
    end
end

function plot_cpp_replay_actuator_differences(t, matlab_results, cpp_results, channel_labels)
    names = common_result_names(matlab_results, cpp_results);
    if isempty(names)
        return;
    end

    colors = lines(numel(names));
    figure('Name', 'SHC09 C++ replay actuator difference', 'Color', 'w');

    for ch = 1:numel(channel_labels)
        subplot(numel(channel_labels), 1, ch);
        hold on;

        for i = 1:numel(names)
            name = names{i};
            n = common_sample_count(t, matlab_results.(name), cpp_results.(name));
            good = common_success_mask(matlab_results.(name), cpp_results.(name), n);
            du = cpp_results.(name).u(ch, 1:n) - matlab_results.(name).u(ch, 1:n);
            du_line = du;
            du_line(~good) = NaN;
            plot(t(1:n), du_line, 'Color', colors(i, :), ...
                'LineWidth', 0.7, 'DisplayName', name);

            if any(~good)
                plot(t(~good), du(~good), 'x', 'Color', colors(i, :), ...
                    'LineWidth', 0.7, 'HandleVisibility', 'off');
            end
        end

        yline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        grid on;
        ylabel(sprintf('du%d %s', ch, channel_labels{ch}), 'Interpreter', 'none');

        if ch == 1
            title('C++ replay 对比: u_{cpp} - u_{matlab}', 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 3, 'Interpreter', 'none');
        end
        if ch == numel(channel_labels)
            xlabel('Time (s)');
        end
    end
end

function plot_cpp_replay_tracking_differences(t, matlab_results, cpp_results, B_eval, axis_labels)
    names = common_result_names(matlab_results, cpp_results);
    if isempty(names)
        return;
    end

    colors = lines(numel(names));
    figure('Name', 'SHC09 C++ replay tracking difference', 'Color', 'w');

    for ax = 1:numel(axis_labels)
        subplot(numel(axis_labels), 1, ax);
        hold on;

        for i = 1:numel(names)
            name = names{i};
            n = common_sample_count(t, matlab_results.(name), cpp_results.(name));
            good = common_success_mask(matlab_results.(name), cpp_results.(name), n);
            du = cpp_results.(name).u(:, 1:n) - matlab_results.(name).u(:, 1:n);
            dBu = B_eval(ax, :) * du;
            dBu_line = dBu;
            dBu_line(~good) = NaN;
            plot(t(1:n), dBu_line, 'Color', colors(i, :), ...
                'LineWidth', 0.7, 'DisplayName', name);

            if any(~good)
                plot(t(~good), dBu(~good), 'x', 'Color', colors(i, :), ...
                    'LineWidth', 0.7, 'HandleVisibility', 'off');
            end
        end

        yline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        grid on;
        ylabel(axis_labels{ax}, 'Interpreter', 'none');

        if ax == 1
            title('C++ replay 对比: B_{par}*(u_{cpp} - u_{matlab})', 'Interpreter', 'none');
            legend('Location', 'best', 'NumColumns', 3, 'Interpreter', 'none');
        end
        if ax == numel(axis_labels)
            xlabel('Time (s)');
        end
    end
end

function names = common_result_names(matlab_results, cpp_results)
    matlab_names = fieldnames(matlab_results);
    cpp_names = fieldnames(cpp_results);
    names = matlab_names(ismember(matlab_names, cpp_names));
end

function n = common_sample_count(t, matlab_result, cpp_result)
    n = min([ ...
        numel(t), ...
        size(matlab_result.u, 2), ...
        size(cpp_result.u, 2), ...
        numel(matlab_result.errout), ...
        numel(cpp_result.errout)]);
end

function good = common_success_mask(matlab_result, cpp_result, n)
    good = matlab_result.errout(1:n) == 0 ...
        & cpp_result.errout(1:n) == 0 ...
        & all(isfinite(matlab_result.u(:, 1:n)), 1) ...
        & all(isfinite(cpp_result.u(:, 1:n)), 1);
end

function report_cpp_plot_skipped(matlab_results, cpp_results)
    if isempty(fieldnames(cpp_results))
        fprintf('C++ replay comparison plots skipped: cpp_results is empty.\n');
        return;
    end

    matlab_names = fieldnames(matlab_results);
    cpp_names = fieldnames(cpp_results);
    fprintf('C++ replay comparison plots skipped: no common method names. MATLAB=[%s], C++=[%s]\n', ...
        strjoin(matlab_names.', ', '), strjoin(cpp_names.', ', '));
end
