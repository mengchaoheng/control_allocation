%% Control allocation benchmark from PX4 ULog topics
%
% 本脚本只做离线控制分配对比，模型统一写成：
%
%     v = B*u
%
% 其中 6 维 v 的行顺序固定为 PX4 ControlAllocator 顺序：
%
%     v(1:3) = 力矩 [Mx; My; Mz]
%     v(4:6) = 力   [Fx; Fy; Fz]
%
% -------------------------------------------------------------------------
% 主流程
% -------------------------------------------------------------------------
% 1. 在本文件开头设置日志、B、算法、C++、绘图选项。
%
% 2. 读取/生成 COMMAND_MAT
%    parse_control_allocation_flight_data.m 负责从原始 PX4 topic 构造：
%      v_sp_t                : vehicle_torque_setpoint_0 timestamp, seconds
%      log_v_sp_instances{i} : 6 x N, [Mx My Mz Fx Fy Fz]
%      u_px4                 : actuator_motors_0 + actuator_servos_0
%
% 3. 按样本直接分配
%      v_sp_i = allocator_instance_v_sp{inst_idx}(:, idx);
%      u_i    = run_allocator(method, v_sp_i(rows), B(rows,:), umin, umax);
%      v_out  = B*u_i;
%
% 4. 保存、打印、画图
%    结果保存到 RESULT_MAT，然后调用 print_control_allocation_comparison.m
%    和 plot_control_allocation_results.m。
%
% -------------------------------------------------------------------------
% 参数说明
% -------------------------------------------------------------------------
% LOG_DIR:
%   放 .ulg 日志的目录。可以是绝对路径或相对路径。
%
% LOG_FILE:
%   .ulg 文件名，例如 "04_08_18.ulg"。
%
% MODEL_NAME:
%   只用于结果标记和打印，例如 "df4"、"shc09"、"shw09"。
%
% PARSE_LOG:
%   true  : 重新从 LOG_PATH 解析 ULog，覆盖 COMMAND_MAT。
%   false : 直接加载已有 COMMAND_MAT。第一次运行或改 parser 后设 true。
%
% SAMPLE_RANGE:
%   []        : 使用全部样本。
%   1:200     : 使用前 200 个样本。
%   500:900   : 使用指定样本段。
%
% METHODS_TO_RUN:
%   {'inv','pca_dir',...} : 只跑列出的算法。
%   ALL_METHODS            : 跑下面列出的全部 MATLAB allocator。
%
% USE_RESTORING:
%   true  : allocator raw u 之后调用 restoring_cpp(B,u,umin,umax)。
%   false : 只看 raw allocator 输出。
%
% NORMALIZE_B:
%   true  : 先把执行器真实幅值并入列 B*diag(u_abs_max)，再用 px4_normalize_B()
%           把 B 转成 PX4 normalized effectiveness。
%   false : 直接用你写的原始 B。
%
% u_abs_max:
%   每个 B case 里的执行器真实最大幅值。allocator 输出 u 仍是归一化值。
%   电机 0..max：u_abs_max = max；舵面左右对称：u_abs_max = abs(max)。
%   若不写，默认全 1。
%
% NORMALIZE_V_SP:
%   true  : 对 v_sp 使用和 B 相同的 D 矩阵缩放。适合物理 B + 物理 v 一起测试。
%   false : 保持日志 v_sp 原值。适合直接回放 PX4 allocator 输入。
%
% NORMALIZE_RPY:
%   true  : px4_normalize_B() 按 PX4 roll/pitch/yaw 归一化逻辑处理前三轴。
%   false : 不对前三轴做 RPY normalize。
%
% RUN_PRINT_AFTER_COMPARE:
%   true/false，是否跑完后调用 print_control_allocation_comparison.m。
%
% RUN_PLOT_AFTER_COMPARE:
%   true/false，是否跑完后调用 plot_control_allocation_results.m。
%
% RUN_CPP_ALLOCATOR:
%   true  : 同时调用 alloc_cpp/build/ca_offline_benchmark，并把 C++ 输出
%           读回到 results(case).alg 里，名称为 cpp_pca_dir/cpp_pca_dpscaled。
%   false : 只跑 MATLAB allocator。
%
% CPP_ALLOCATOR_TO_PLOT:
%   {}                                      : 不画 C++ 输出。
%   {'cpp_pca_dir'}                         : 只画 C++ pca_dir。
%   {'cpp_pca_dpscaled'}                    : 只画 C++ pca_dpscaled。
%   {'cpp_pca_dir','cpp_pca_dpscaled'}      : 两个 C++ 输出都画。
%
% CPP_ALLOCATOR_EXE:
%   C++ 离线分配程序路径。当前程序输入/输出格式是：
%     results/cpp_input/inst_0_B.csv, inst_0_Y.csv, ...
%     results/cpp_output/pca_dir_u.csv, pca_dpscaled_u.csv, ...
%
% PLOT_SAVE_PNG / PLOT_SAVE_FIG / PLOT_SHOW_FIGURES:
%   控制 plot 脚本保存 png、保存 fig、是否显示 MATLAB 图窗。
%
% PLOT_NORMALIZED:
%   true  : 图里显示归一化的 u、v_sp、B*u、residual。
%   false : 图里 u 乘回 u_abs_max；B*u 用 D 还原。
%           v_sp 只有在 NORMALIZE_V_SP=true 时才用 D 还原。

clc; clear; close all;

tool_dir = fileparts(mfilename('fullpath'));
addpath(genpath(tool_dir));

%% 1. User settings
LOG_DIR = "/Users/mch/Proj/Mac_DF/PX4-Autopilot/build/px4_sitl_default/rootfs/log/2026-06-04";
LOG_FILE = "04_08_18.ulg";
LOG_PATH = fullfile(LOG_DIR, LOG_FILE);
MODEL_NAME = "df4";

RESULT_DIR = fullfile(tool_dir, 'results');
[~, log_stem] = fileparts(char(LOG_PATH));
COMMAND_MAT = fullfile(RESULT_DIR, [log_stem '_command_data.mat']);
RESULT_MAT = fullfile(RESULT_DIR, [log_stem '_allocation_compare_results.mat']);
FIGURE_DIR = fullfile(RESULT_DIR, [log_stem '_figures']);

PARSE_LOG = false;
SAMPLE_RANGE = [];

ALL_METHODS = { ...
    'inv', ...
    'pca_dir', ...
    'pca_dpscaled', ...
    'pca_prio', ...
    'lib_lpwrap_db', ...
    'lib_lpwrap_dbinf', ...
    'lib_lpwrap_dir', ...
    'lib_lpwrap_dpscaled', ...
    'lib_lpwrap_mo', ...
    'lib_lpwrap_sb', ...
    'lib_lpwrap_par_db', ...
    'lib_lpwrap_par_dbinf', ...
    'lib_lpwrap_par_dir', ...
    'lib_lpwrap_par_dpscaled', ...
    'lib_lpwrap_par_mo', ...
    'lib_lpwrap_par_sb', ...
    'lib_lpwrap_incre', ...
    'lib_cgiwrap', ...
    'lib_dawrap', ...
    'lib_vjawrap', ...
    'wls', ...
    'wls_gen', ...
    'dir_linprog', ...
    'dir_linprog_re', ...
    'dir_linprog_re_bound', ...
    'use_lp_lib', ...
    'allocator_dir_lpwrap_4'};
METHODS_TO_RUN = {'inv', 'pca_dir', 'pca_dpscaled', 'wls'};
% METHODS_TO_RUN = ALL_METHODS;
USE_RESTORING = true;
NORMALIZE_B = true;
NORMALIZE_V_SP = false;
NORMALIZE_RPY = true;

RUN_PRINT_AFTER_COMPARE = true;
RUN_PLOT_AFTER_COMPARE = true;
PLOT_SAVE_PNG = true;
PLOT_SAVE_FIG = false;
PLOT_SHOW_FIGURES = usejava('desktop');
PLOT_NORMALIZED = true;

RUN_CPP_ALLOCATOR = true;
CPP_ALLOCATOR_TO_PLOT = {'cpp_pca_dir', 'cpp_pca_dpscaled'};
CPP_ALLOCATOR_EXE = fullfile(tool_dir, 'alloc_cpp', 'build', 'ca_offline_benchmark');

B_CASES = {};
% B_CASES{end+1} = make_df4_single_case();
% B_CASES{end+1} = make_df4_split_case();
B_CASES{end+1} = make_shc09_single_case();
B_CASES{end+1} = make_shc09_split_case();

%% 2. Load raw ULog topic tables
if PARSE_LOG || ~isfile(COMMAND_MAT)
    run(fullfile(tool_dir, 'parse_control_allocation_flight_data.m'));
end

load(COMMAND_MAT);

%% 3. Select samples and run allocation
if isempty(SAMPLE_RANGE)
    sample_idx = 1:numel(v_sp_t);
else
    sample_idx = SAMPLE_RANGE;
end

v_sp_t = v_sp_t(sample_idx);
for i = 1:numel(log_v_sp_instances)
    log_v_sp_instances{i} = log_v_sp_instances{i}(:, sample_idx);
end
u_px4 = u_px4(sample_idx, :);

flightData = struct();
flightData.model = char(MODEL_NAME);
flightData.t = v_sp_t(:);
flightData.v_sp = log_v_sp_instances{1}';
flightData.u_px4 = u_px4;
flightData.u_px4_source = u_px4_source;

benchmark_tic = tic;
results = repmat(empty_case_result(), 1, numel(B_CASES));
methods_to_run = cellfun(@(x) lower(char(x)), cellstr(METHODS_TO_RUN), 'UniformOutput', false);
METHODS_TO_PLOT = [methods_to_run, cellstr(CPP_ALLOCATOR_TO_PLOT)];

fprintf('\nControl allocation benchmark\n');
fprintf('  setpoint mat: %s\n', COMMAND_MAT);
fprintf('  methods     : %s\n', strjoin(methods_to_run, ', '));
fprintf('  plot        : %s\n', strjoin(METHODS_TO_PLOT, ', '));
fprintf('  restoring   : %s\n', yesno(USE_RESTORING));
fprintf('  normalize B : %s\n', yesno(NORMALIZE_B));
fprintf('  normalize v : %s\n', yesno(NORMALIZE_V_SP));
fprintf('  C++ bridge  : %s\n', yesno(RUN_CPP_ALLOCATOR));

for case_idx = 1:numel(B_CASES)
    cfg = B_CASES{case_idx};
    allocator_instances = prepare_allocator_instances(cfg, NORMALIZE_B, NORMALIZE_RPY);
    allocator_instance_v_sp = map_log_v_sp_to_allocator_instances(log_v_sp_instances, numel(allocator_instances));
    allocator_instance_v_sp = normalize_v_sp_if_needed(allocator_instance_v_sp, allocator_instances, NORMALIZE_V_SP);
    v_sp = zeros(size(allocator_instance_v_sp{1}));
    for inst_idx = 1:numel(allocator_instance_v_sp)
        v_sp = v_sp + allocator_instance_v_sp{inst_idx};
    end

    fprintf('\n[Case %d/%d] %s\n', case_idx, numel(B_CASES), cfg.name);
    fprintf('  B: 6 x %d, instances: %d, samples: %d\n', ...
        sum(arrayfun(@(x) size(x.B, 2), allocator_instances)), numel(allocator_instances), size(v_sp, 2));

    alg = repmat(empty_algorithm_result(), 1, numel(methods_to_run));

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};
        alg(method_idx) = run_method_on_log(method, allocator_instances, allocator_instance_v_sp, v_sp, USE_RESTORING);

        fprintf('  %-22s total=%7.4fs  avg=%8.2fus  rms(v-Bu)=%9.4g  max=%9.4g  fail=%d\n', ...
            alg(method_idx).name, alg(method_idx).elapsed_s, ...
            alg(method_idx).avg_us_per_sample, ...
            alg(method_idx).rmse_residual, alg(method_idx).max_abs_residual, ...
            alg(method_idx).fail_count);
    end

    cpp_process_wall_s = nan;
    if RUN_CPP_ALLOCATOR
        [cpp_alg, cpp_process_wall_s] = run_cpp_allocator_case(CPP_ALLOCATOR_EXE, RESULT_DIR, case_idx, cfg.name, ...
            allocator_instances, allocator_instance_v_sp, v_sp, USE_RESTORING, flightData);
        alg = [alg cpp_alg]; %#ok<AGROW>
    end

    [alg, u_inv_ref, inv_reference_name] = attach_inv_reference(alg);

    result = empty_case_result();
    result.name = cfg.name;
    result.B = join_B(allocator_instances);
    result.D = join_D(allocator_instances);
    result.umin = vertcat(allocator_instances.umin);
    result.umax = vertcat(allocator_instances.umax);
    result.u_abs_max = vertcat(allocator_instances.u_abs_max);
    result.rows = find(any(abs(result.B) > 1e-9, 2))';
    result.cpp_process_wall_s = cpp_process_wall_s;
    result.t = v_sp_t(:);
    result.v_sp = v_sp';
    result.u_inv_ref = u_inv_ref;
    result.inv_reference_name = inv_reference_name;
    result.alg = compare_with_px4(alg, flightData);
    results(case_idx) = result;
end

benchmark_elapsed_s = toc(benchmark_tic);

meta = struct();
meta.MODEL_NAME = char(MODEL_NAME);
meta.LOG_PATH = LOG_PATH;
meta.COMMAND_MAT = COMMAND_MAT;
meta.SAMPLE_RANGE = SAMPLE_RANGE;
meta.METHODS_TO_RUN = METHODS_TO_RUN;
meta.CPP_ALLOCATOR_TO_PLOT = CPP_ALLOCATOR_TO_PLOT;
meta.METHODS_TO_PLOT = METHODS_TO_PLOT;
meta.USE_RESTORING = USE_RESTORING;
meta.NORMALIZE_B = NORMALIZE_B;
meta.NORMALIZE_V_SP = NORMALIZE_V_SP;
meta.PLOT_NORMALIZED = PLOT_NORMALIZED;
meta.RUN_CPP_ALLOCATOR = RUN_CPP_ALLOCATOR;
meta.CPP_ALLOCATOR_EXE = CPP_ALLOCATOR_EXE;
meta.benchmark_elapsed_s = benchmark_elapsed_s;

save(RESULT_MAT, 'results', 'flightData', 'meta', 'COMMAND_MAT', 'LOG_PATH', '-v7.3');

fprintf('\nSaved results:\n  %s\n', RESULT_MAT);
fprintf('Benchmark wall time: %.4f s\n', benchmark_elapsed_s);

if RUN_PRINT_AFTER_COMPARE
    run(fullfile(tool_dir, 'print_control_allocation_comparison.m'));
end

if RUN_PLOT_AFTER_COMPARE
    SAVE_PNG = PLOT_SAVE_PNG;
    SAVE_FIG = PLOT_SAVE_FIG;
    SHOW_FIGURES = PLOT_SHOW_FIGURES;
    run(fullfile(tool_dir, 'plot_control_allocation_results.m'));
end

%% B cases
function cfg = make_df4_single_case()
B_force = [0; 0; 0; 0; 0; -27.8];
B_moment = make_ducted_fan_moment_B("cross4", [0.01149 0.01153 0.00487], 0.167, 0.069, 1.0);
servo_max = 40 * pi / 180;

cfg.name = 'df4_single';
cfg.B = [B_force, [B_moment; zeros(3, 4)]];
cfg.umin = [0; -ones(4, 1)];
cfg.umax = [1;  ones(4, 1)];
cfg.u_abs_max = [1; servo_max * ones(4, 1)];
end

function cfg = make_df4_split_case()
B_force = [0; 0; 0; 0; 0; -27.8];
B_moment = make_ducted_fan_moment_B("cross4", [0.01149 0.01153 0.00487], 0.167, 0.069, 1.0);
servo_max = 40 * pi / 180;

cfg.name = 'df4_split';
cfg.B = {B_force, [B_moment; zeros(3, 4)]};
cfg.umin = {0, -ones(4, 1)};
cfg.umax = {1,  ones(4, 1)};
cfg.u_abs_max = {1, servo_max * ones(4, 1)};
end

function cfg = make_shc09_single_case()
B_force = [0; 0; 0; 0; 0; -6.5];
B_moment = make_ducted_fan_moment_B("hex6", [0.0438 0.0436 0.005006], 0.267, 0.066, 1.0);
servo_max = 40 * pi / 180;

cfg.name = 'shc09_single';
cfg.B = [B_force, [B_moment; zeros(3, 6)]];
cfg.umin = [0; -ones(6, 1)];
cfg.umax = [1;  ones(6, 1)];
cfg.u_abs_max = [1; servo_max * ones(6, 1)];
end

function cfg = make_shc09_split_case()
B_force = [0; 0; 0; 0; 0; -6.5];
B_moment = make_ducted_fan_moment_B("hex6", [0.0438 0.0436 0.005006], 0.267, 0.066, 1.0);
servo_max = 40 * pi / 180;

cfg.name = 'shc09_split';
cfg.B = {B_force, [B_moment; zeros(3, 6)]};
cfg.umin = {0, -ones(6, 1)};
cfg.umax = {1,  ones(6, 1)};
cfg.u_abs_max = {1, servo_max * ones(6, 1)};
end

function B_moment = make_ducted_fan_moment_B(layout, inertia_diag, l1, l2, k_surface)
I = diag(double(inertia_diag(:)));
d = 60 * pi / 180;

switch lower(string(layout))
    case {"cross4", "df4"}
        P = [-l1,  0,   l1,  0;
              0,  -l1,  0,   l1;
              l2,  l2,  l2,  l2];

    case {"hex6", "df6", "shc09"}
        P = [-l1, -cos(d)*l1,  cos(d)*l1,  l1,  cos(d)*l1, -cos(d)*l1;
              0,   sin(d)*l1,  sin(d)*l1,  0,  -sin(d)*l1, -sin(d)*l1;
              l2,  l2,         l2,         l2,  l2,         l2];

    otherwise
        error('Unknown B layout: %s', layout);
end

B_moment = I \ P * k_surface;
end

%% B and instance preparation
function allocator_instances = prepare_allocator_instances(cfg, normalize_B, normalize_rpy)
if iscell(cfg.B)
    B_cell = cfg.B;
    umin_cell = cfg.umin;
    umax_cell = cfg.umax;
else
    B_cell = {cfg.B};
    umin_cell = {cfg.umin};
    umax_cell = {cfg.umax};
end

if isfield(cfg, 'u_abs_max')
    if iscell(cfg.u_abs_max)
        u_abs_max_cell = cfg.u_abs_max;
    else
        u_abs_max_cell = {cfg.u_abs_max};
    end
else
    u_abs_max_cell = cell(size(B_cell));
end

allocator_instances = repmat(struct('B', [], 'B_alloc', [], 'D', [], 'rows', [], ...
    'umin', [], 'umax', [], 'u_abs_max', [], 'cols', []), 1, numel(B_cell));
next_col = 1;

for i = 1:numel(B_cell)
    B_raw = force_6_rows(B_cell{i});
    u_abs_max = resolve_u_abs_max(u_abs_max_cell{i}, size(B_raw, 2));
    B_raw = B_raw * diag(u_abs_max);

    if normalize_B
        [D, B, ~, ~] = px4_normalize_B(B_raw, normalize_rpy);
        B = force_6_rows(B);
    else
        D = eye(6);
        B = B_raw;
    end

    rows = find(any(abs(B) > 1e-9, 2))';
    cols = next_col:(next_col + size(B, 2) - 1);
    next_col = next_col + size(B, 2);

    allocator_instances(i).B = B;
    allocator_instances(i).B_alloc = B(rows, :);
    allocator_instances(i).D = D;
    allocator_instances(i).rows = rows;
    allocator_instances(i).umin = double(umin_cell{i}(:));
    allocator_instances(i).umax = double(umax_cell{i}(:));
    allocator_instances(i).u_abs_max = u_abs_max;
    allocator_instances(i).cols = cols;
end
end

function u_abs_max = resolve_u_abs_max(u_abs_max, actuator_count)
if isempty(u_abs_max)
    u_abs_max = ones(actuator_count, 1);
else
    u_abs_max = double(u_abs_max(:));
end

if isscalar(u_abs_max)
    u_abs_max = repmat(u_abs_max, actuator_count, 1);
end

if numel(u_abs_max) ~= actuator_count
    error('u_abs_max must be scalar or have %d entries.', actuator_count);
end
end

function allocator_instance_v_sp = map_log_v_sp_to_allocator_instances(log_v_sp_instances, num_allocator_instances)
% 日志 instance 数和当前 B instance 数一致：逐个喂给对应 B。
if num_allocator_instances == numel(log_v_sp_instances)
    allocator_instance_v_sp = log_v_sp_instances;
    return;
end

% 日志 2 个 instance，当前 B 1 个 instance：
% 合成 v，前三行力矩来自日志 instance 1，后三行力来自日志 instance 0。
if num_allocator_instances == 1 && numel(log_v_sp_instances) >= 2
    v_sp = zeros(size(log_v_sp_instances{1}));
    v_sp(1:3, :) = log_v_sp_instances{2}(1:3, :);
    v_sp(4:6, :) = log_v_sp_instances{1}(4:6, :);
    allocator_instance_v_sp = {v_sp};
    return;
end

% 日志 1 个 instance，当前 B 2 个 instance：
% 默认第 0 个 B 是力分配，吃 rows 4:6；第 1 个 B 是力矩分配，吃 rows 1:3。
if num_allocator_instances == 2 && numel(log_v_sp_instances) == 1
    v_sp = log_v_sp_instances{1};
    allocator_instance_v_sp = {zeros(size(v_sp)), zeros(size(v_sp))};
    allocator_instance_v_sp{1}(4:6, :) = v_sp(4:6, :);
    allocator_instance_v_sp{2}(1:3, :) = v_sp(1:3, :);
    return;
end

error('Cannot map %d log allocation instance(s) to %d B allocation instance(s).', ...
    numel(log_v_sp_instances), num_allocator_instances);
end

function allocator_instance_v_sp = normalize_v_sp_if_needed(allocator_instance_v_sp, allocator_instances, normalize_v_sp)
if ~normalize_v_sp
    return;
end

for i = 1:numel(allocator_instance_v_sp)
    allocator_instance_v_sp{i} = allocator_instances(i).D * allocator_instance_v_sp{i};
end
end

%% Allocation
function alg = run_method_on_log(method, allocator_instances, allocator_instance_v_sp, v_sp, use_restoring)
N = size(v_sp, 2);
m_total = sum(arrayfun(@(x) size(x.B, 2), allocator_instances));
u = nan(N, m_total);
v_alloc = nan(N, 6);
residual = nan(N, 6);
restore_s = 0;
fail_count = 0;

elapsed_tic = tic;

for idx = 1:N
    v_alloc_sample = zeros(6, 1);
    sample_failed = false;

    for inst_idx = 1:numel(allocator_instances)
        inst = allocator_instances(inst_idx);
        rows = inst.rows;
        v_sp_i = allocator_instance_v_sp{inst_idx}(:, idx);

        try
            u_i = run_allocator(method, v_sp_i(rows), inst.B_alloc, inst.umin, inst.umax);
            [u_i, restore_dt] = apply_restoring_if_needed(inst.B_alloc, u_i, inst.umin, inst.umax, use_restoring);
            v_alloc_sample = v_alloc_sample + inst.B * u_i;
        catch
            u_i = nan(size(inst.B, 2), 1);
            restore_dt = 0;
            sample_failed = true;
            v_alloc_sample(:) = nan;
        end

        u(idx, inst.cols) = u_i';
        restore_s = restore_s + restore_dt;
    end

    v_alloc(idx, :) = v_alloc_sample';
    residual(idx, :) = v_sp(:, idx)' - v_alloc_sample';
    fail_count = fail_count + double(sample_failed);
end

elapsed_s = toc(elapsed_tic);
[rms_residual, max_residual] = finite_rmse_max(residual);

alg = empty_algorithm_result();
alg.name = char(method);
alg.method = char(method);
alg.restoring_mode = ternary(use_restoring, 'restored', 'raw');
alg.u = u;
alg.v_achieved = v_alloc;
alg.residual = residual;
alg.elapsed_s = elapsed_s;
alg.allocator_s = elapsed_s - restore_s;
alg.restore_s = restore_s;
alg.avg_us_per_sample = 1e6 * alg.allocator_s / max(N, 1);
alg.avg_restore_us_per_sample = 1e6 * restore_s / max(N, 1);
alg.rmse_residual = rms_residual;
alg.max_abs_residual = max_residual;
alg.fail_count = fail_count;
alg.fallback_count = 0;
end

function u = run_allocator(method, v, B, umin, umax)
[k, m] = size(B);
global NumU
NumU = m;

switch lower(char(method))
    case 'inv'
        u = pinv(B) * v;

    case 'pca_dir'
        [~, u] = evalc('DP_LPCA(v, B, umin, umax, 100);');

    case 'pca_dpscaled'
        [~, u] = evalc('DPscaled_LPCA(v, B, umin, umax, 100);');

    case 'pca_prio'
        [~, u] = evalc('DP_LPCA_prio(zeros(k, 1), v, B, umin, umax, 100);');

    case 'lib_lpwrap_db'
        [~, u] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin, umax, 0));');

    case 'lib_lpwrap_dbinf'
        [~, u] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin, umax, 1));');

    case 'lib_lpwrap_dir'
        [~, u] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin, umax, 2));');

    case 'lib_lpwrap_dpscaled'
        [~, u] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin, umax, 3));');

    case 'lib_lpwrap_mo'
        [~, u] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin, umax, 4));');

    case 'lib_lpwrap_sb'
        [~, u] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin, umax, 5));');

    case 'lib_lpwrap_par_db'
        [~, u] = evalc('LPwrap_par(make_lpwrap_in_mat(B, v, umin, umax, 0), v, m);');

    case 'lib_lpwrap_par_dbinf'
        [~, u] = evalc('LPwrap_par(make_lpwrap_in_mat(B, v, umin, umax, 1), v, m);');

    case 'lib_lpwrap_par_dir'
        [~, u] = evalc('LPwrap_par(make_lpwrap_in_mat(B, v, umin, umax, 2), v, m);');

    case 'lib_lpwrap_par_dpscaled'
        [~, u] = evalc('LPwrap_par(make_lpwrap_in_mat(B, v, umin, umax, 3), v, m);');

    case 'lib_lpwrap_par_mo'
        [~, u] = evalc('LPwrap_par(make_lpwrap_in_mat(B, v, umin, umax, 4), v, m);');

    case 'lib_lpwrap_par_sb'
        [~, u] = evalc('LPwrap_par(make_lpwrap_in_mat(B, v, umin, umax, 5), v, m);');

    case 'lib_lpwrap_incre'
        u0 = zeros(m, 1);
        [~, du] = evalc('LPwrap(make_lpwrap_in_mat(B, v, umin - u0, umax - u0, 2));');
        u = du(:) + u0;

    case 'lib_cgiwrap'
        [~, u] = evalc('CGIwrap(make_book_wrapper_in_mat(B, v, umin, umax));');

    case 'lib_dawrap'
        [~, u] = evalc('DAwrap(make_book_wrapper_in_mat(B, v, umin, umax));');

    case 'lib_vjawrap'
        [~, u] = evalc('VJAwrap(make_book_wrapper_in_mat(B, v, umin, umax));');

    case 'wls'
        Wv = eye(k);
        Wu = eye(m);
        ud = zeros(m, 1);
        gam = 1e6;
        u0 = (umin + umax) / 2;
        W0 = zeros(m, 1);
        [~, u] = evalc('wls_alloc(B, v, umin, umax, Wv, Wu, ud, gam, u0, W0, 100);');

    case 'wls_gen'
        Wv = eye(k);
        Wu = eye(m);
        ud = zeros(m, 1);
        gam = 1e6;
        u0 = zeros(m, 1);
        W0 = zeros(m, 1);
        [~, u] = evalc('wls_alloc_gen(B, v, umin, umax, Wv, Wu, ud, gam, u0, W0, 100, m);');

    case 'dir_linprog'
        [~, u, ~] = evalc('dir_alloc_linprog(B, v, umin, umax, 1e4);');

    case 'dir_linprog_re'
        [~, u, ~] = evalc('dir_alloc_linprog_re(B, v, umin, umax);');

    case 'dir_linprog_re_bound'
        [~, u, ~] = evalc('dir_alloc_linprog_re_bound(B, v, umin, umax, 1e4);');

    case 'use_lp_lib'
        [~, u, ~] = evalc('use_LP_lib(B, v, umin, umax);');

    case 'allocator_dir_lpwrap_4'
        if m ~= 4
            error('allocator_dir_LPwrap_4 only supports m=4, current m=%d.', m);
        end
        [~, u, ~, ~] = evalc('allocator_dir_LPwrap_4(single(B), single(v), single(umin), single(umax));');
        u = double(u(:));

    otherwise
        error('Unknown allocator method: %s', method);
end

u = clamp_u(u, umin, umax);
end

function [u, elapsed_s] = apply_restoring_if_needed(B, u_raw, umin, umax, use_restoring)
if use_restoring
    restore_tic = tic;
    u = restoring_cpp(B, u_raw(:), umin, umax);
    elapsed_s = toc(restore_tic);
else
    u = u_raw(:);
    elapsed_s = 0;
end
end

%% C++ offline allocator bridge
function [cpp_alg, cpp_wall_s] = run_cpp_allocator_case(cpp_exe, result_dir, ~, ~, ...
    allocator_instances, allocator_instance_v_sp, v_sp, use_restoring, flightData)
cpp_alg = repmat(empty_algorithm_result(), 1, 2);
cpp_wall_s = nan;

if ~isfile(cpp_exe)
    warning('C++ allocator executable not found: %s', cpp_exe);
    cpp_alg = struct([]);
    return;
end

cpp_input_dir = fullfile(result_dir, 'cpp_input');
cpp_output_dir = fullfile(result_dir, 'cpp_output');

if ~isfolder(cpp_input_dir)
    mkdir(cpp_input_dir);
end
if ~isfolder(cpp_output_dir)
    mkdir(cpp_output_dir);
end

N = size(v_sp, 2);
total_u_dim = sum(arrayfun(@(x) size(x.B, 2), allocator_instances));
writematrix([N, total_u_dim, numel(allocator_instances), double(use_restoring)], ...
    fullfile(cpp_input_dir, 'case_meta.csv'));

for inst_idx = 1:numel(allocator_instances)
    inst = allocator_instances(inst_idx);
    prefix = fullfile(cpp_input_dir, sprintf('inst_%d', inst_idx - 1));
    writematrix(inst.B, [prefix '_B.csv']);
    writematrix(allocator_instance_v_sp{inst_idx}', [prefix '_Y.csv']);
    writematrix(inst.umin(:)', [prefix '_umin.csv']);
    writematrix(inst.umax(:)', [prefix '_umax.csv']);
end

cmd = sprintf('"%s" "%s" "%s"', cpp_exe, cpp_input_dir, cpp_output_dir);
cpp_wall_tic = tic;
[status, cmdout] = system(cmd);
cpp_wall_s = toc(cpp_wall_tic);

if status ~= 0
    warning('C++ allocator failed:\n%s', cmdout);
    cpp_alg = struct([]);
    return;
end

fprintf('\nC++ allocator bridge:\n');
fprintf('  input : %s\n', cpp_input_dir);
fprintf('  output: %s\n', cpp_output_dir);
fprintf('  wall  : %.4fs\n', cpp_wall_s);

cpp_alg(1) = read_cpp_algorithm_result(cpp_output_dir, 'pca_dir', N, flightData);
cpp_alg(2) = read_cpp_algorithm_result(cpp_output_dir, 'pca_dpscaled', N, flightData);

for i = 1:numel(cpp_alg)
    fprintf('  %-22s total=%7.4fs  avg=%8.2fus  rms(v-Bu)=%9.4g  max=%9.4g  fail=%d  fb=%d\n', ...
        cpp_alg(i).name, cpp_alg(i).elapsed_s, cpp_alg(i).avg_us_per_sample, ...
        cpp_alg(i).rmse_residual, cpp_alg(i).max_abs_residual, ...
        cpp_alg(i).fail_count, cpp_alg(i).fallback_count);
end
end

function alg = read_cpp_algorithm_result(cpp_output_dir, method, N, flightData)
u = fixed_rows(readmatrix(fullfile(cpp_output_dir, [method '_u.csv'])), N);
v_achieved = fixed_rows(readmatrix(fullfile(cpp_output_dir, [method '_v_achieved.csv'])), N);
residual = fixed_rows(readmatrix(fullfile(cpp_output_dir, [method '_residual.csv'])), N);
timing = readmatrix(fullfile(cpp_output_dir, [method '_timing.csv']));

[rms_residual, max_residual] = finite_rmse_max(residual);

alg = empty_algorithm_result();
alg.name = ['cpp_' method];
alg.method = ['cpp_' method];
alg.restoring_mode = 'restored';
alg.u = u;
alg.v_achieved = v_achieved;
alg.residual = residual;
alg.elapsed_s = timing(1);
alg.allocator_s = timing(2);
alg.restore_s = timing(3);
alg.avg_us_per_sample = 1e6 * timing(2) / max(N, 1);
alg.avg_restore_us_per_sample = 1e6 * timing(3) / max(N, 1);
alg.rmse_residual = rms_residual;
alg.max_abs_residual = max_residual;
alg.fail_count = round(timing(4));
alg.fallback_count = round(timing(5));

if isfield(flightData, 'u_px4') && size(flightData.u_px4, 2) == size(u, 2)
    [alg.rmse_vs_px4, alg.max_abs_vs_px4] = finite_rmse_max(u - flightData.u_px4(1:size(u, 1), :));
end
end

function X = fixed_rows(X, N)
if size(X, 1) < N
    X(end+1:N, :) = nan;
elseif size(X, 1) > N
    X = X(1:N, :);
end
end

function IN_MAT = make_lpwrap_in_mat(B, v, umin, umax, lp_method)
[~, m] = size(B);
IN_MAT = [B v(:);
          umin(:)' 0;
          umax(:)' 0;
          ones(1, m) lp_method];
end

function IN_MAT = make_book_wrapper_in_mat(B, v, umin, umax)
IN_MAT = make_lpwrap_in_mat(B, v, umin, umax, 0);
end

%% Result helpers
function [alg, u_inv_ref, inv_reference_name] = attach_inv_reference(alg)
u_inv_ref = [];
inv_reference_name = "";
idx = find(strcmp({alg.method}, 'inv'), 1);

if isempty(idx)
    return;
end

u_inv_ref = alg(idx).u;
inv_reference_name = string(alg(idx).name);

for i = 1:numel(alg)
    [alg(i).rmse_vs_inv, alg(i).max_abs_vs_inv] = finite_rmse_max(alg(i).u - u_inv_ref);
end
end

function alg = compare_with_px4(alg, flightData)
for i = 1:numel(alg)
    alg(i).rmse_vs_px4 = nan;
    alg(i).max_abs_vs_px4 = nan;

    if isfield(flightData, 'u_px4') && size(flightData.u_px4, 2) == size(alg(i).u, 2)
        N = min(size(flightData.u_px4, 1), size(alg(i).u, 1));
        [alg(i).rmse_vs_px4, alg(i).max_abs_vs_px4] = ...
            finite_rmse_max(alg(i).u(1:N, :) - flightData.u_px4(1:N, :));
    end
end
end

function result = empty_case_result()
result = struct('name', '', 'B', [], 'D', [], 'umin', [], 'umax', [], 'rows', [], ...
    'u_abs_max', [], 'cpp_process_wall_s', nan, 't', [], 'v_sp', [], ...
    'u_inv_ref', [], 'inv_reference_name', "", 'alg', empty_algorithm_result());
end

function alg = empty_algorithm_result()
alg = struct('name', '', 'method', '', 'restoring_mode', '', ...
    'u', [], 'v_achieved', [], 'residual', [], ...
    'elapsed_s', nan, 'allocator_s', nan, 'restore_s', nan, ...
    'avg_us_per_sample', nan, 'avg_restore_us_per_sample', nan, ...
    'rmse_residual', nan, 'max_abs_residual', nan, ...
    'rmse_vs_inv', nan, 'max_abs_vs_inv', nan, ...
    'rmse_vs_px4', nan, 'max_abs_vs_px4', nan, ...
    'fail_count', 0, 'fallback_count', 0);
end

function B = join_B(allocator_instances)
B = zeros(6, sum(arrayfun(@(x) size(x.B, 2), allocator_instances)));

for i = 1:numel(allocator_instances)
    B(:, allocator_instances(i).cols) = allocator_instances(i).B;
end
end

function D = join_D(allocator_instances)
D = eye(6);

for i = 1:numel(allocator_instances)
    rows = allocator_instances(i).rows;
    D(rows, rows) = allocator_instances(i).D(rows, rows);
end
end

function X = force_6_rows(X)
X = double(X);
X(~isfinite(X)) = 0;

if size(X, 1) < 6
    X(end+1:6, :) = 0;
elseif size(X, 1) > 6
    error('Expected [Mx My Mz Fx Fy Fz] rows, got %dx%d.', size(X, 1), size(X, 2));
end
end

function u = clamp_u(u, umin, umax)
u = min(max(u(:), umin(:)), umax(:));
end

function [rmse_value, max_value] = finite_rmse_max(X)
v = X(isfinite(X));

if isempty(v)
    rmse_value = nan;
    max_value = nan;
else
    rmse_value = sqrt(mean(v.^2));
    max_value = max(abs(v));
end
end

function text = yesno(tf)
text = ternary(tf, 'yes', 'no');
end

function value = ternary(condition, a, b)
if condition
    value = a;
else
    value = b;
end
end
