function plot_all_results(t, u_reference, results, channel_labels, B_eval, y_eval, axis_labels, reference_name, u_reference_label)
% 打开控制分配测试的所有关键图，不保存 PNG。
%   1) 每个执行器通道 u_reference 和各算法 u 的输出曲线。
%   2) 每个执行器通道 u_algorithm-u_reference 的误差曲线。
%   3) 每个控制轴 B*u-y 的 tracking error。
%   4) 每个执行器通道相对 reference_name 的输出差异。
    if nargin < 9 || isempty(u_reference_label)
        u_reference_label = 'u_reference';
    end

    plot_actuator_outputs(t, u_reference, results, channel_labels, u_reference_label);
    plot_actuator_reference_line_errors(t, u_reference, results, channel_labels, u_reference_label);
    plot_tracking_errors(t, results, B_eval, y_eval, axis_labels);
    plot_actuator_reference_errors(t, results, channel_labels, reference_name);
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

function plot_tracking_errors(t, results, B_eval, y_eval, axis_labels)
    names = fieldnames(results);
    figure('Name', 'SHC09 tracking error B*u-y', 'Color', 'w');

    for ax = 1:numel(axis_labels)
        subplot(numel(axis_labels), 1, ax);
        hold on;

        for i = 1:numel(names)
            u = results.(names{i}).u;
            err = B_eval(ax, :) * u - y_eval(ax, :);
            plot(t, err, 'LineWidth', 0.7, 'DisplayName', names{i});
        end

        yline(0, 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off');
        grid on;
        ylabel(axis_labels{ax}, 'Interpreter', 'none');

        if ax == 1
            title('控制量跟踪误差: B_{par}*u - y_{par}', 'Interpreter', 'none');
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
