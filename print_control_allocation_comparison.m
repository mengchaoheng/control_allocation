%% Print minimal control-allocation benchmark result
% 只打印三类信息：
%   1. 每个 B case 的基本信息；
%   2. 每个算法的时间、residual、fail/fallback 计数；
%   3. 算法之间的 u 和 B*u 差异。
%
% 主参考是每个当前 B case 离线重算的 inv。PX4 日志 actuator 输出只有在列数
% 和当前 B 完全一致时才打印，避免把不同机型/不同 actuator 定义拿来硬比。

%% 读取结果
% 这个文件既可以被 main_control_allocation_benchmark.m 用 run() 调用，
% 也可以单独运行。
%   - 如果 main 已经在工作区里放好了 results/flightData，就直接打印；
%   - 如果单独运行，就像 plot 文件一样自动加载最新的结果 .mat。
tool_dir = fileparts(mfilename('fullpath'));

if exist('results', 'var') ~= 1 || exist('flightData', 'var') ~= 1
    if ~exist('RESULT_MAT', 'var')
        RESULT_MAT = "";
    else
        RESULT_MAT = string(RESULT_MAT);
    end

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

    fprintf('\nLoaded result file:\n  %s\n', RESULT_MAT);
end

axis_names = ["Mx", "My", "Mz", "Fx", "Fy", "Fz"];

fprintf('\n============================================================\n');
fprintf('Control allocation benchmark\n');
fprintf('Model   : %s\n', flightData.model);
fprintf('Samples : %d\n', size(flightData.v_sp, 1));

if exist('meta', 'var') == 1 && isfield(meta, 'benchmark_elapsed_s')
    fprintf('Wall    : %.4f s\n', meta.benchmark_elapsed_s);
end

fprintf('============================================================\n');

%% 0. PX4 online allocator runtime from allocation_value
print_online_allocation_runtime(flightData);

%% 1. 每个 B / 每个算法
fprintf('\n[Per B / per algorithm]\n');

for case_idx = 1:numel(results)
    r = results(case_idx);
    active_axes = axis_names(r.rows);

    fprintf('\nCase %d: %s\n', case_idx, r.name);
    fprintf('  B           : %dx%d\n', size(r.B, 1), size(r.B, 2));
    fprintf('  active axes : %s\n', join_or_none(active_axes));
    fprintf('  samples     : %d\n', size(r.v_sp, 1));
    fprintf('  inv ref     : %s\n', string_or_none(get_field_or(r, 'inv_reference_name', "")));

    if isfield(r, 'cpp_process_wall_s') && isfinite(r.cpp_process_wall_s)
        fprintf('  C++ wall    : %.4f s, process start + CSV I/O + C++ loops\n', r.cpp_process_wall_s);
    end

    fprintf('\n');
    fprintf('  method              total(s)  avg_alloc(us)  restore(us)  rms(v-Bu)   max(v-Bu)   fail  fb\n');
    fprintf('  ------------------  --------  -------------  -----------  ----------  ----------  ----  ---\n');

    for alg_idx = 1:numel(r.alg)
        a = r.alg(alg_idx);

        fprintf('  %-18s  %8.4f  %13.2f  %11.2f  %10.4g  %10.4g  %4d  %3d\n', ...
            a.name, a.elapsed_s, a.avg_us_per_sample, ...
            1e6 * a.restore_s / max(size(r.v_sp, 1), 1), ...
            a.rmse_residual, a.max_abs_residual, a.fail_count, a.fallback_count);
    end

    fprintf('\n');
    fprintf('  Offline inv reference on the same B:\n');
    fprintf('  method              rms(u-u_inv)  max(u-u_inv)\n');
    fprintf('  ------------------  ------------  ------------\n');

    for alg_idx = 1:numel(r.alg)
        a = r.alg(alg_idx);

        fprintf('  %-18s  %12.4g  %12.4g\n', ...
            a.name, a.rmse_vs_inv, a.max_abs_vs_inv);
    end

    ref_name = reference_output_name(flightData);
    has_ref = any(arrayfun(@(a) isfinite(a.rmse_vs_px4), r.alg));

    if has_ref
        fprintf('\n');
        fprintf('  PX4 actuator output (%s):\n', ref_name);
        fprintf('  method              rms(u-u_px4)  max(u-u_px4)\n');
        fprintf('  ------------------  ------------  ------------\n');

        for alg_idx = 1:numel(r.alg)
            a = r.alg(alg_idx);

            if isfinite(a.rmse_vs_px4)
                fprintf('  %-18s  %12.4g  %12.4g\n', ...
                    a.name, a.rmse_vs_px4, a.max_abs_vs_px4);
            end
        end
    else
        ref_dim = reference_output_dim(flightData);

        if ref_dim > 0
            fprintf('\n');
            fprintf('  PX4 actuator output skipped: %s has %d column(s), current B has %d actuator(s).\n', ...
                ref_name, ref_dim, size(r.B, 2));
        end
    end
end

%% 2. 同一个 B，不同算法
fprintf('\n[Same B, different algorithms]\n');

for case_idx = 1:numel(results)
    r = results(case_idx);

    if numel(r.alg) < 2
        continue;
    end

    fprintf('\nCase: %s\n', r.name);
    fprintf('  method A            method B            rms(du)     max(du)     rms(dBu)    max(dBu)\n');
    fprintf('  ------------------  ------------------  ----------  ----------  ----------  ----------\n');

    for i = 1:(numel(r.alg) - 1)
        for j = (i + 1):numel(r.alg)
            a = r.alg(i);
            b = r.alg(j);
            [u_rmse, u_max] = finite_rmse_max(a.u - b.u);
            [v_rmse, v_max] = finite_rmse_max(a.v_achieved - b.v_achieved);

            fprintf('  %-18s  %-18s  %10.4g  %10.4g  %10.4g  %10.4g\n', ...
                a.name, b.name, u_rmse, u_max, v_rmse, v_max);
        end
    end
end

%% 3. 同一个算法，不同 B
fprintf('\n[Same algorithm, different B]\n');

if numel(results) >= 2
    method_names = string({results(1).alg.name});

    for method_idx = 1:numel(method_names)
        fprintf('\nAlgorithm: %s\n', method_names(method_idx));
        fprintf('  case A              case B              rms(dBu)    max(dBu)\n');
        fprintf('  ------------------  ------------------  ----------  ----------\n');

        for i = 1:(numel(results) - 1)
            for j = (i + 1):numel(results)
                a = results(i).alg(method_idx);
                b = results(j).alg(method_idx);
                N = min(size(a.v_achieved, 1), size(b.v_achieved, 1));
                [v_rmse, v_max] = finite_rmse_max(a.v_achieved(1:N, :) - b.v_achieved(1:N, :));

                fprintf('  %-18s  %-18s  %10.4g  %10.4g\n', ...
                    results(i).name, results(j).name, v_rmse, v_max);
            end
        end
    end
else
    fprintf('  only one B case enabled\n');
end

fprintf('\nDone.\n');

function text = join_or_none(values)
if isempty(values)
    text = 'none';
else
    text = char(join(values, ' '));
end
end

function [rmse_v, max_v] = finite_rmse_max(x)
v = x(isfinite(x));

if isempty(v)
    rmse_v = nan;
    max_v = nan;
else
    rmse_v = sqrt(mean(v.^2));
    max_v = max(abs(v));
end
end

function value = get_field_or(s, name, default_value)
if isfield(s, name)
    value = s.(name);
else
    value = default_value;
end
end

function text = string_or_none(value)
text = string(value);

if strlength(text) == 0
    text = "none";
end

text = char(text);
end

function dim = reference_output_dim(flightData)
dim = 0;

if ~isempty(flightData.u_px4)
    dim = size(flightData.u_px4, 2);
end
end

function name = reference_output_name(flightData)
if isfield(flightData, 'u_px4_source') && strlength(string(flightData.u_px4_source)) > 0
    name = char(string(flightData.u_px4_source));
else
    name = 'u_px4';
end
end

function print_online_allocation_runtime(flightData)
fprintf('\n[PX4 online allocation runtime]\n');

if ~isfield(flightData, 'allocation_runtime') || isempty(flightData.allocation_runtime)
    fprintf('  no allocation_value timing in result file\n');
    return;
end

runtime = flightData.allocation_runtime;
fprintf('  instance  method          samples  avg(us)  prep/core/post(us)\n');
fprintf('  --------  --------------  -------  -------  ------------------\n');

for i = 1:numel(runtime)
    s = runtime(i);
    fprintf('  %8d  %-14s  %7d  %7.2f  %6.2f/%6.2f/%6.2f\n', ...
        s.instance, s.method_name, s.samples, s.mean_us, ...
        s.prepare_mean_us, s.core_mean_us, s.post_mean_us);
end
end

function result_mat = select_latest_result_mat(result_dir)
% 自动选择最新结果文件，方便单独运行 print 脚本。
files = dir(fullfile(result_dir, '*_allocation_compare_results.mat'));

if isempty(files)
    error('目录里没有 *_allocation_compare_results.mat: %s', result_dir);
end

[~, newest_idx] = max([files.datenum]);
result_mat = string(fullfile(files(newest_idx).folder, files(newest_idx).name));
end
