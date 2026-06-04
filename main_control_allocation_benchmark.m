%% Minimal control-allocation benchmark
% 目标很窄：用真实飞行日志里的 allocator 输入 y(t)，测试用户手写 B 下的不同分配算法。
%
% 主流程只有四步：
%   1. 读日志里 ControlAllocator::Run() 的输入：
%        y = [vehicle_torque_setpoint.xyz, vehicle_thrust_setpoint.xyz]
%   2. 用户在本文件里手写一个或多个 B 实例，每个 B 都整理成 6 x n：
%        行顺序 [Mx My Mz Fx Fy Fz]，列顺序就是执行器顺序。
%   3. 预处理 B：
%        可选 px4_normalize_B() 单位化；
%        找非零行 rows；
%        之后算法只取 y(rows) 和 B(rows,:)。
%   4. 对每个 B case、每个算法、每个样本、每个实例分配，最后拼成一个 u。

clc;clear all; close all;

%% 1. 路径和日志
tool_dir = fileparts(mfilename('fullpath'));
external_allocator_repo = string(tool_dir);
addpath(tool_dir);

if isfolder(external_allocator_repo)
    % 只加入本脚本实际调用的算法目录，避免 genpath 扫大量无关 demo/LP 库。
    addpath(char(external_allocator_repo));  % px4_normalize_B.m
    addpath(fullfile(char(external_allocator_repo), 'PCA'));
    addpath(fullfile(char(external_allocator_repo), 'restoring'));
    addpath(fullfile(char(external_allocator_repo), 'control_allocation_lib', 'qcat', 'QCAT', 'qcat'));
end

% 路径只保留两个入口：
%   LOG_DIR    : 放 .ulg 的文件夹；
%   RESULT_DIR : 放解析后的 command MAT、分配结果 MAT、绘图文件。
FORCE_REPARSE = false;   % true = 不管 MAT 是否存在，都重新从 .ulg 解析。

LOG_DIR = "/Users/mch/Proj/Mac_DF/PX4-Autopilot/build/px4_sitl_default/rootfs/log/2026-06-04";
LOG_FILE = "04_08_18.ulg";
LOG_PATH = fullfile(LOG_DIR, LOG_FILE);
MODEL_NAME = "df4";
COMMAND_INSTANCE = 0;

RESULT_DIR = fullfile(tool_dir, 'results');
[~, log_stem] = fileparts(char(LOG_PATH));
COMMAND_MAT = fullfile(RESULT_DIR, [log_stem '_command_data.mat']);
RESULT_MAT = fullfile(RESULT_DIR, [log_stem '_allocation_compare_results.mat']);
FIGURE_DIR = fullfile(RESULT_DIR, [log_stem '_figures']);

% 二选一：WINDOW_MODE="sample" 使用 TEST_SAMPLE_WINDOW；
%        WINDOW_MODE="time"   使用 TEST_TIME_WINDOW_S；
%        WINDOW_MODE="all"    跑全部样本。
WINDOW_MODE = "all";
TEST_SAMPLE_WINDOW = 1:200;
TEST_TIME_WINDOW_S = [];

%% 2. 用户手写 B
% 加新机型/新 B，只复制一个 case_cfg block。
%
% 单实例：
%   case_cfg.B = B;
%   case_cfg.umin = umin;
%   case_cfg.umax = umax;
%
% 多实例：
%   case_cfg.B = {B0, B1};
%   case_cfg.umin = {umin0, umin1};
%   case_cfg.umax = {umax0, umax1};
%
% 每个实例都必须是 6 x n。某个实例不控制的轴直接填 0。
B_CASES = {};

% df4 split 示例：
%   instance 0: motor，主要产生 Fz；
%   instance 1: 四个舵面，产生 Mx/My/Mz。
% 力矩部分用 make_ducted_fan_moment_B() 构造，避免以后每个机型都手抄符号矩阵。
B0 = [
      0
      0
      0
      0
      0
    -27.8
];

df4_I = [0.01149 0.01153 0.00487];
df4_l1 = 0.167;
df4_l2 = 0.069;
df4_k_moment = 1.0;
B1_moment = make_ducted_fan_moment_B("cross4", df4_I, df4_l1, df4_l2, df4_k_moment);
B1 = [B1_moment; zeros(3, size(B1_moment, 2))];

case_cfg = struct();
case_cfg.name = "df4_split";
case_cfg.B = {B0, B1};
case_cfg.umin = {0, -ones(4, 1)};
case_cfg.umax = {1,  ones(4, 1)};
case_cfg.normalize_rpy = [true true];
B_CASES{end+1} = case_cfg;

% 单实例 B 示例：
case_cfg = struct();
case_cfg.name = "df4_single";
case_cfg.B = [B0 B1];
case_cfg.umin = [0; -ones(4, 1)];
case_cfg.umax = [1;  ones(4, 1)];
case_cfg.normalize_rpy = true;
B_CASES{end+1} = case_cfg;

% B case 选择：'all' = 全部；数字 = B_CASES 序号；字符串/cell = case_cfg.name。
B_CASE_SELECTION = 'all';

% 常用力矩 B 构造模板：
%   df6:
%     Bm = make_ducted_fan_moment_B("hex6", [0.01149 0.01153 0.00487], ...
%          0.167, 0.069, 1.0, 60*pi/180);
%
%   SHC09:
%     Bm = make_ducted_fan_moment_B("hex6", [0.0438 0.0436 0.005006], ...
%          0.267, 0.066, 1.0, 60*pi/180);
%
%   SHW09_vtol hover/MC 8 舵面：
%     wing_yaw_torque = k_wing * L_3;   % 若 MC 阶段翼面不参与，可填 [] 或 0。
%     Bm = make_ducted_fan_moment_B("hex6_wing2", [I_x I_y I_z], ...
%          L_1, L_2, k_cs, 60*pi/180, wing_yaw_torque);
%
% 得到 Bm 后，放入某个 allocation instance：
%     B_surface = [Bm; zeros(3, size(Bm, 2))];  % 只控制 Mx/My/Mz，不控制 Fx/Fy/Fz。

%% 3. 算法和预处理开关
METHOD_CATALOG = {'inv', 'pca_dir', 'pca_dpscaled', 'pca_prio', 'wls'};
METHOD_SELECTION = {'inv', 'pca_dir', 'pca_dpscaled', 'wls'};
RESTORING_SELECTION = "restored";    % "raw" / "restored" / "both"
SIMPLEX_BACKEND = "original";         % "original" 更快；"tiebreak" 用于和当前 C++ guarded tie-break 对齐。
WARMUP_SAMPLE_COUNT = 100;            % MATLAB JIT warmup，不计入算法耗时；设 0 可看冷启动时间。
FALLBACK_METHOD = "inv";
ENABLE_FALLBACK = true;
ENFORCE_U_LIMITS = true;

% 单位化只处理 B，不改日志输入 y。
USE_UNIT_ALLOCATION = true;
DEFAULT_NORMALIZE_RPY = true;
ZERO_ROW_TOL = 1e-9;

% main 运行完可以自动打印/画图，不需要手动再运行 print/plot。
RUN_PRINT_AFTER_COMPARE = true;
RUN_PLOT_AFTER_COMPARE = true;
PLOT_SAVE_PNG = true;
PLOT_SAVE_FIG = false;
PLOT_SHOW_FIGURES = usejava('desktop');  % 桌面 MATLAB 弹图；batch/headless 只保存。

%% 4. 解析或读取日志输入
% 第一次运行时 RESULT_DIR 里没有 command MAT，main 会自动解析日志；
% 以后默认直接加载 MAT，避免把 ulog2csv 的时间混进分配测试。
if FORCE_REPARSE || ~isfile(COMMAND_MAT)
    if ~isfolder(RESULT_DIR)
        mkdir(RESULT_DIR);
    end

    run(fullfile(tool_dir, 'parse_control_allocation_flight_data.m'));
end

if ~isfile(COMMAND_MAT)
    error('解析后仍找不到 COMMAND_MAT: %s。请检查 LOG_PATH 或 ulog2csv。', COMMAND_MAT);
end

S = load(COMMAND_MAT, 'flightData');
flightData = select_window(S.flightData, WINDOW_MODE, TEST_SAMPLE_WINDOW, TEST_TIME_WINDOW_S);
log_Y = read_log_inputs(flightData);

%% 5. 预处理 B，然后进入分配循环
if ~isfolder(RESULT_DIR)
    mkdir(RESULT_DIR);
end

methods_to_run = resolve_selection(METHOD_SELECTION, METHOD_CATALOG, 'allocator method');
B_CASES = select_B_cases(B_CASES, B_CASE_SELECTION);
jobs = make_jobs(methods_to_run, RESTORING_SELECTION);
RESTORING_AVAILABLE = exist('restoring_cpp', 'file') == 2;
TIE_OPTS = make_simplex_options(SIMPLEX_BACKEND);
results_cell = cell(1, numel(B_CASES));
benchmark_tic = tic;

fprintf('\nControl allocation offline benchmark\n');
fprintf('  command mat : %s\n', COMMAND_MAT);
fprintf('  log input   : %d instance(s), %d sample(s)\n', numel(log_Y), size(log_Y{1}, 1));
fprintf('  methods     : %s\n', strjoin(string(methods_to_run), ', '));
fprintf('  restoring   : %s\n', RESTORING_SELECTION);
fprintf('  simplex     : %s\n', SIMPLEX_BACKEND);
fprintf('  warmup      : %d sample(s), not counted in per-method timing\n', WARMUP_SAMPLE_COUNT);

for case_idx = 1:numel(B_CASES)
    B_case = prepare_B_case(B_CASES{case_idx}, USE_UNIT_ALLOCATION, ...
        DEFAULT_NORMALIZE_RPY, ZERO_ROW_TOL);

    % 输入映射规则：
    %   日志实例数 == B 实例数：一一对应。
    %   日志单实例、B 多实例：instance0 吃力/F，instance1 吃力矩/M。
    %   日志多实例、B 单实例：所有日志实例相加后给这个大 B。
    Y_cell = make_instance_inputs(log_Y, B_case.inst);
    N = min(cellfun(@(Y) size(Y, 1), Y_cell));
    y_command = sum_instance_inputs(Y_cell, N);

    print_case_header(case_idx, numel(B_CASES), B_case, flightData, N);

    alg_cell = cell(1, numel(jobs));
    total_u_dim = numel(B_case.umin);

    for job_idx = 1:numel(jobs)
        job = jobs(job_idx);
        check_method_available(job.method);
        warmup_allocator_job(job, B_case.inst, Y_cell, WARMUP_SAMPLE_COUNT, ...
            ENFORCE_U_LIMITS, ENABLE_FALLBACK, FALLBACK_METHOD, TIE_OPTS, RESTORING_AVAILABLE);

        u_all = nan(N, total_u_dim);
        y_achieved = nan(N, 6);
        residual = nan(N, 6);

        primary_time = 0;
        fallback_time = 0;
        restore_time = 0;
        fail_count = 0;
        fallback_count = 0;
        alg_tic = tic;

        % 核心分配循环：样本 -> B 实例。
        for sample_idx = 1:N
            y_out = zeros(6, 1);
            sample_failed = false;

            for inst_idx = 1:numel(B_case.inst)
                inst = B_case.inst(inst_idx);
                y_full = Y_cell{inst_idx}(sample_idx, :)';

                % 只有当前 B 实例能产生的控制轴才进入算法。
                rows = inst.rows;
                if isempty(rows)
                    u_raw = zeros(inst.m, 1);
                    u = u_raw;
                    B_eval = inst.B;
                    did_fallback = false;
                else
                    [B_run, y_run] = allocator_input(job.method, inst, y_full);
                    umin = inst.umin;
                    umax = inst.umax;

                    if ~ENFORCE_U_LIMITS
                        umin = -inf(size(umin));
                        umax = inf(size(umax));
                    end

                    try
                        [u_raw, did_fallback, primary_dt, fallback_dt] = ...
                            run_allocator_with_fallback(job.method, y_run, B_run, ...
                            y_full, inst.B, umin, umax, ENABLE_FALLBACK, FALLBACK_METHOD, TIE_OPTS);

                        primary_time = primary_time + primary_dt;
                        fallback_time = fallback_time + fallback_dt;

                        [u, restore_dt] = apply_optional_restoring(inst.B_eff, u_raw, ...
                            umin, umax, job.use_restoring, RESTORING_AVAILABLE);
                        restore_time = restore_time + restore_dt;

                        B_eval = inst.B;
                    catch ME
                        u_raw = nan(inst.m, 1);
                        u = u_raw;
                        B_eval = inst.B;
                        did_fallback = false;
                        sample_failed = true;
                    end
                end

                u_all(sample_idx, inst.cols) = u(:)';
                y_out = y_out + B_eval * u;
                fallback_count = fallback_count + double(did_fallback);

                if ~all(isfinite(u)) || any(u < inst.umin - 1e-7) || any(u > inst.umax + 1e-7)
                    sample_failed = true;
                end
            end

            y_achieved(sample_idx, :) = y_out';
            residual(sample_idx, :) = y_command(sample_idx, :) - y_out';
            fail_count = fail_count + double(sample_failed);
        end

        allocator_time = primary_time + fallback_time;
        elapsed_s = toc(alg_tic);
        residual_values = residual(isfinite(residual));
        [rmse_vs_px4, max_abs_vs_px4] = compare_with_px4_output(u_all, flightData);
        alg_cell{job_idx} = struct( ...
            'name', job.name, ...
            'method', job.method, ...
            'restoring_mode', job.restoring_mode, ...
            'u', u_all, ...
            'y_achieved', y_achieved, ...
            'residual', residual, ...
            'elapsed_s', elapsed_s, ...
            'allocator_s', allocator_time, ...
            'restore_s', restore_time, ...
            'avg_us_per_sample', 1e6 * allocator_time / max(N, 1), ...
            'avg_restore_us_per_sample', 1e6 * restore_time / max(N, 1), ...
            'rmse_residual', sqrt(mean(residual_values.^2)), ...
            'max_abs_residual', max(abs(residual_values)), ...
            'rmse_vs_px4', rmse_vs_px4, ...
            'max_abs_vs_px4', max_abs_vs_px4, ...
            'fail_count', fail_count, ...
            'fallback_count', fallback_count);

        fprintf('  %-14s total=%6.4fs avg_alloc=%7.2fus residual_rms=%9.4g fb=%d fail=%d\n', ...
            job.name, elapsed_s, alg_cell{job_idx}.avg_us_per_sample, ...
            alg_cell{job_idx}.rmse_residual, fallback_count, fail_count);
    end

    results_cell{case_idx} = struct( ...
        'name', B_case.name, ...
        'B', B_case.B, ...
        'umin', B_case.umin, ...
        'umax', B_case.umax, ...
        'rows', B_case.rows, ...
        't', time_vector(flightData, N), ...
        'control_sp', y_command, ...
        'alg', [alg_cell{:}]);
end

results = [results_cell{:}];
benchmark_elapsed_s = toc(benchmark_tic);

meta = struct();
meta.MODEL_NAME = flightData.model;
meta.COMMAND_MAT = COMMAND_MAT;
meta.LOG_PATH = LOG_PATH;
meta.FORCE_REPARSE = FORCE_REPARSE;
meta.WINDOW_MODE = WINDOW_MODE;
meta.TEST_SAMPLE_WINDOW = TEST_SAMPLE_WINDOW;
meta.TEST_TIME_WINDOW_S = TEST_TIME_WINDOW_S;
meta.METHODS_TO_RUN = methods_to_run;
meta.RESTORING_SELECTION = RESTORING_SELECTION;
meta.RESTORING_AVAILABLE = RESTORING_AVAILABLE;
meta.SIMPLEX_BACKEND = SIMPLEX_BACKEND;
meta.WARMUP_SAMPLE_COUNT = WARMUP_SAMPLE_COUNT;
meta.USE_UNIT_ALLOCATION = USE_UNIT_ALLOCATION;
meta.DEFAULT_NORMALIZE_RPY = DEFAULT_NORMALIZE_RPY;
meta.benchmark_elapsed_s = benchmark_elapsed_s;
meta.px4_reference_note = 'u_px4 is actuator_motors+actuator_servos interpolated onto command samples.';

save(RESULT_MAT, 'results', 'flightData', 'meta', 'COMMAND_MAT', 'LOG_PATH', '-v7.3');

fprintf('\nSaved results:\n  %s\n', RESULT_MAT);
fprintf('Benchmark wall time: %.4f s\n', benchmark_elapsed_s);

if RUN_PRINT_AFTER_COMPARE
    run(fullfile(tool_dir, 'print_control_allocation_comparison.m'));
end

if RUN_PLOT_AFTER_COMPARE
    SAVE_PNG = PLOT_SAVE_PNG; %#ok<NASGU>
    SAVE_FIG = PLOT_SAVE_FIG; %#ok<NASGU>
    SHOW_FIGURES = PLOT_SHOW_FIGURES; %#ok<NASGU>
    run(fullfile(tool_dir, 'plot_control_allocation_results.m'));
end

%% B preprocessing
function B_case = prepare_B_case(cfg, use_unit_allocation, default_normalize_rpy, zero_tol)
% 把用户输入的单实例/多实例 B 统一成 6 x n cell，并为每个实例保存：
%   B_raw : 原始 B；
%   B     : 单位化后的 B；
%   rows  : 非零控制轴索引；
%   B_eff : B(rows,:)；分配算法实际吃这个矩阵。
if iscell(cfg.B)
    B_cell = cfg.B;
    umin_cell = cfg.umin;
    umax_cell = cfg.umax;
else
    B_cell = {cfg.B};
    umin_cell = {cfg.umin};
    umax_cell = {cfg.umax};
end

num_inst = numel(B_cell);
inst = repmat(struct('B_raw', [], 'B', [], 'B_eff', [], ...
    'rows', [], 'umin', [], 'umax', [], 'cols', [], 'm', 0), 1, num_inst);

B_all = [];
umin_all = [];
umax_all = [];
next_col = 1;

for i = 1:num_inst
    B_raw = force_6_rows(B_cell{i});
    umin = double(umin_cell{i}(:));
    umax = double(umax_cell{i}(:));

    if size(B_raw, 2) ~= numel(umin) || size(B_raw, 2) ~= numel(umax)
        error('%s instance %d: B 列数必须等于 umin/umax 长度。', cfg.name, i - 1);
    end

    normalize_rpy = case_normalize_rpy(cfg, i, default_normalize_rpy);

    if use_unit_allocation
        B = unitize_B(B_raw, normalize_rpy);
    else
        B = B_raw;
    end

    rows = find(any(abs(B) > zero_tol, 2));
    cols = next_col:(next_col + size(B, 2) - 1);
    next_col = next_col + size(B, 2);

    inst(i).B_raw = B_raw;
    inst(i).B = B;
    inst(i).B_eff = B(rows, :);
    inst(i).rows = rows(:)';
    inst(i).umin = umin;
    inst(i).umax = umax;
    inst(i).cols = cols;
    inst(i).m = size(B, 2);

    B_all = [B_all B]; %#ok<AGROW>
    umin_all = [umin_all; umin]; %#ok<AGROW>
    umax_all = [umax_all; umax]; %#ok<AGROW>
end

B_case = struct();
B_case.name = char(cfg.name);
B_case.inst = inst;
B_case.B = B_all;
B_case.umin = umin_all;
B_case.umax = umax_all;
B_case.rows = unique([inst.rows]);
B_case.unit_allocation_enabled = use_unit_allocation;
end

function B_moment = make_ducted_fan_moment_B(layout, inertia_diag, l1, l2, k_surface, d_rad, wing_yaw_torque)
% 构造这类 ducted-fan/tailsitter 机型最常用的力矩效应矩阵。
%
% 输出：
%   B_moment 是 3 x n，只包含 [Mx; My; Mz] 三个力矩轴。
%   main 里再用 [B_moment; zeros(3,n)] 放进 6 x n allocation instance。
%
% 输入：
%   layout          "cross4"      : df4 四个尾舵；
%                   "hex6"       : df6 / SHC09 六个环形尾舵；
%                   "hex6_wing2" : SHW09_vtol 六个环形尾舵 + 两个翼面 yaw。
%   inertia_diag    [I_x I_y I_z]。
%   l1              roll/pitch 力臂。
%   l2              yaw 力臂。
%   k_surface       舵面偏转到侧向力的等效增益，例如 dF_d_delta_h_cs。
%                   若你已经在测试里使用简化无量纲 B，可直接填 1。
%   d_rad           六舵布局夹角，默认 60 deg；cross4 不使用它。
%   wing_yaw_torque SHW09 翼面 yaw 预惯量力矩增益 k_wing * L_3。
%                   填 [] 或 0 表示两个翼面暂不参与。
if nargin < 6 || isempty(d_rad)
    d_rad = 60 * pi / 180;
end

if nargin < 7
    wing_yaw_torque = [];
end

I = diag(double(inertia_diag(:)));
layout = lower(string(layout));

switch layout
    case {"cross4", "df4"}
        P = [-l1,  0,   l1,  0;
               0, -l1,  0,   l1;
              l2,  l2,  l2,  l2];
        B_moment = I \ P * k_surface;

    case {"hex6", "df6", "shc09"}
        P = [-l1, -cos(d_rad)*l1,  cos(d_rad)*l1,  l1,  cos(d_rad)*l1, -cos(d_rad)*l1;
               0,  sin(d_rad)*l1,  sin(d_rad)*l1,  0,  -sin(d_rad)*l1, -sin(d_rad)*l1;
              l2,  l2,             l2,             l2,  l2,             l2];
        B_moment = I \ P * k_surface;

    case {"hex6_wing2", "shw09_vtol"}
        P = [-l1, -cos(d_rad)*l1,  cos(d_rad)*l1,  l1,  cos(d_rad)*l1, -cos(d_rad)*l1;
               0,  sin(d_rad)*l1,  sin(d_rad)*l1,  0,  -sin(d_rad)*l1, -sin(d_rad)*l1;
              l2,  l2,             l2,             l2,  l2,             l2];
        B_moment = [I \ P * k_surface, zeros(3, 2)];

        if ~isempty(wing_yaw_torque) && abs(wing_yaw_torque) > 0
            % SHW09_vtol 里两个翼面主要按相反方向贡献 yaw。
            B_moment(3, 7:8) = [wing_yaw_torque / inertia_diag(3), ...
                               -wing_yaw_torque / inertia_diag(3)];
        end

    otherwise
        error('Unknown ducted fan moment layout: %s', layout);
end
end

function normalize_rpy = case_normalize_rpy(cfg, idx, default_normalize_rpy)
% normalize_rpy 可以对整个 case 写一个标量，也可以按实例写数组。
if isfield(cfg, 'normalize_rpy') && ~isempty(cfg.normalize_rpy)
    raw = cfg.normalize_rpy;
else
    raw = default_normalize_rpy;
end

if numel(raw) >= idx
    normalize_rpy = logical(raw(idx));
else
    normalize_rpy = logical(raw(end));
end
end

function B_unit = unitize_B(B_raw, normalize_rpy)
% 单位化固定使用现有 px4_normalize_B.m，不在这里维护第二套公式。
if exist('px4_normalize_B', 'file') ~= 2
    error('px4_normalize_B.m not found. Check external_allocator_repo at top of main.');
end

[~, B_unit] = px4_normalize_B(B_raw, normalize_rpy);
B_unit = force_6_rows(B_unit);
end

%% Log input mapping
function log_Y = read_log_inputs(flightData)
% parser 已经把 vehicle_torque_setpoint + vehicle_thrust_setpoint 整理成 control_sp。
if isfield(flightData, 'setpoint_instances') && ~isempty(flightData.setpoint_instances)
    log_Y = cell(1, numel(flightData.setpoint_instances));

    for i = 1:numel(log_Y)
        log_Y{i} = force_6_cols(flightData.setpoint_instances(i).control_sp);
    end
else
    log_Y = {force_6_cols(flightData.control_sp)};
end
end

function Y_cell = make_instance_inputs(log_Y, inst)
% 输入实例数不一致时只做当前需要的简单转换：
%   日志单实例 -> B 多实例：force 给第一个 B，torque 给第二个 B。
%   日志多实例 -> B 单实例：所有日志实例相加给一个大 B。
num_log = numel(log_Y);
num_B = numel(inst);

if num_log == num_B
    Y_cell = log_Y;
    return;
end

if num_log > 1 && num_B == 1
    N = min(cellfun(@(Y) size(Y, 1), log_Y));
    Y_cell = {sum_instance_inputs(log_Y, N)};
    return;
end

if num_log == 1 && num_B > 1
    Y = log_Y{1};
    Y_cell = cell(1, num_B);

    for i = 1:num_B
        Y_cell{i} = zeros(size(Y));
    end

    Y_cell{1}(:, 4:6) = Y(:, 4:6);  % force/thrust -> matrix0
    Y_cell{2}(:, 1:3) = Y(:, 1:3);  % torque/moment -> matrix1
    return;
end

error('不支持日志实例数 %d -> B 实例数 %d 的输入转换。', num_log, num_B);
end

function Y = sum_instance_inputs(Y_cell, N)
% 把多个 6 轴输入相加，得到本次 benchmark 的总 command y。
Y = zeros(N, 6);

for i = 1:numel(Y_cell)
    Y = Y + Y_cell{i}(1:N, :);
end
end

%% Allocation
function jobs = make_jobs(methods, restoring)
% 支持 raw/restored/both，方便对比恢复算法前后的输出。
if ischar(methods) || isstring(methods)
    methods = cellstr(methods);
end

restoring = lower(string(restoring));
if restoring == "both"
    restore_modes = ["raw", "restored"];
else
    restore_modes = restoring;
end

jobs = repmat(struct('method', '', 'name', '', 'use_restoring', false, 'restoring_mode', ''), 1, numel(methods) * numel(restore_modes));
k = 0;

for i = 1:numel(methods)
    for j = 1:numel(restore_modes)
        k = k + 1;
        jobs(k).method = lower(char(methods{i}));
        jobs(k).use_restoring = restore_modes(j) == "restored";
        jobs(k).restoring_mode = char(restore_modes(j));

        if numel(restore_modes) == 1
            jobs(k).name = jobs(k).method;
        else
            jobs(k).name = sprintf('%s_%s', jobs(k).method, restore_modes(j));
        end
    end
end
end

function selected = resolve_selection(selection, catalog, label)
% 和旧 test.m 一样：支持 'all'、数字索引、方法名/case 名。
catalog = cellstr(catalog);

if ischar(selection) || isstring(selection)
    if isscalar(string(selection)) && strcmpi(string(selection), "all")
        selected = catalog;
    else
        selected = cellstr(selection);
    end
elseif isnumeric(selection)
    selected = catalog(selection);
elseif iscell(selection)
    if numel(selection) == 1 && (ischar(selection{1}) || isstring(selection{1})) ...
            && strcmpi(string(selection{1}), "all")
        selected = catalog;
    else
        selected = cellfun(@char, selection, 'UniformOutput', false);
    end
else
    error('Unsupported %s selection type.', label);
end

unknown = setdiff(selected, catalog, 'stable');
if ~isempty(unknown)
    error('Unknown %s: %s. Available: %s', label, strjoin(unknown, ', '), strjoin(catalog, ', '));
end
end

function selected_cases = select_B_cases(B_cases, selection)
names = cellfun(@(c) char(c.name), B_cases, 'UniformOutput', false);
selected_names = resolve_selection(selection, names, 'B case');
selected_cases = cell(1, numel(selected_names));

for i = 1:numel(selected_names)
    idx = find(strcmp(names, selected_names{i}), 1);
    selected_cases{i} = B_cases{idx};
end
end

function tie_opts = make_simplex_options(simplex_backend)
% DP_LPCA/DPscaled_LPCA 的 simplex 后端。
% original 速度更快；tiebreak 用于和 guarded tie-break 版本对照。
switch lower(string(simplex_backend))
    case "original"
        tie_opts = struct('simplex_backend', 'original');
    case "tiebreak"
        tie_opts = struct();
    otherwise
        error('未知 SIMPLEX_BACKEND=%s。只能是 "original" 或 "tiebreak"。', simplex_backend);
end
end

function warmup_allocator_job(job, inst_list, Y_cell, warmup_count, ...
    enforce_limits, enable_fallback, fallback_method, tie_opts, restoring_available)
% MATLAB 第一次跑 DP_LPCA/DPscaled_LPCA 会触发 JIT/函数加载，时间会明显偏大。
% warmup 不保存结果，也不计入 allocator 时间，只用于得到稳定运行时的平均耗时。
if warmup_count <= 0
    return;
end

N = min([warmup_count, cellfun(@(Y) size(Y, 1), Y_cell)]);

for sample_idx = 1:N
    for inst_idx = 1:numel(inst_list)
        inst = inst_list(inst_idx);

        if isempty(inst.rows)
            continue;
        end

        y_full = Y_cell{inst_idx}(sample_idx, :)';
        [B_run, y_run] = allocator_input(job.method, inst, y_full);
        umin = inst.umin;
        umax = inst.umax;

        if ~enforce_limits
            umin = -inf(size(umin));
            umax = inf(size(umax));
        end

        try
            [u_raw, ~, ~, ~] = run_allocator_with_fallback(job.method, ...
                y_run, B_run, y_full, inst.B, umin, umax, ...
                enable_fallback, fallback_method, tie_opts);

            apply_optional_restoring(inst.B_eff, u_raw, umin, umax, job.use_restoring, restoring_available);
        catch
            % 真正的错误会在正式循环里暴露；warmup 不打断 benchmark。
        end
    end
end
end

function check_method_available(method)
% 缺算法文件时直接报错，不在每个样本里重复判断。
switch lower(method)
    case 'pca_dir'
        required = 'DP_LPCA';
    case 'pca_dpscaled'
        required = 'DPscaled_LPCA';
    case 'pca_prio'
        required = 'DP_LPCA_prio';
    case 'wls'
        required = 'wls_alloc';
    otherwise
        required = '';
end

if ~isempty(required) && exist(required, 'file') ~= 2
    error('%s.m not found for method %s.', required, method);
end
end

function [B_run, y_run] = allocator_input(method, inst, y_full)
% inv/wls 使用完整 6 轴 B；PCA/LPCA 只使用有效行 B_eff。
if method_uses_effective_rows(method)
    B_run = inst.B_eff;
    y_run = y_full(inst.rows);
else
    B_run = inst.B;
    y_run = y_full;
end
end

function tf = method_uses_effective_rows(method)
tf = any(strcmpi(method, {'pca_dir', 'pca_dpscaled', 'pca_prio'}));
end

function [u, elapsed_s] = apply_optional_restoring(B, u_raw, umin, umax, use_restoring, restoring_available)
% 是否 restoring 只由一个布尔开关决定；B 满行秩由用户输入时保证。
if use_restoring && restoring_available && ~isempty(B)
    restore_tic = tic;
    u = restoring_cpp(B, u_raw(:), umin, umax);
    elapsed_s = toc(restore_tic);
else
    u = u_raw(:);
    elapsed_s = 0;
end
end

function [u, did_fallback, primary_dt, fallback_dt] = ...
    run_allocator_with_fallback(method, y_run, B_run, ...
    y_full, B_full, umin, umax, enable_fallback, fallback_method, tie_opts)
% PCA/LPCA 在退化有效 B 上没有冗余自由度，直接回退到 inv。
primary_dt = 0;
fallback_dt = 0;
did_fallback = false;

fallback_reason = direct_fallback_reason(method, B_run);

if strlength(fallback_reason) > 0
    if ~enable_fallback
        error('%s cannot run on this B: %s', method, fallback_reason);
    end

    fallback_tic = tic;
    u = run_one_allocator(char(fallback_method), y_full, B_full, umin, umax, tie_opts);
    fallback_dt = toc(fallback_tic);
    did_fallback = true;
    return;
end

try
    primary_tic = tic;
    u = run_one_allocator(method, y_run, B_run, umin, umax, tie_opts);
    primary_dt = toc(primary_tic);
catch ME
    primary_dt = toc(primary_tic);

    if ~enable_fallback
        rethrow(ME);
    end

    fallback_tic = tic;
    u = run_one_allocator(char(fallback_method), y_full, B_full, umin, umax, tie_opts);
    fallback_dt = toc(fallback_tic);
    did_fallback = true;
end
end

function reason = direct_fallback_reason(method, B)
% 有效行维度 k<=1，或执行器数 m<=k 时，没有 PCA/LPCA 可优化的零空间。
reason = "";

if ~method_uses_effective_rows(method)
    return;
end

[k, m] = size(B);

if k <= 1
    reason = sprintf('effective row dimension k=%d <= 1', k);
elseif m <= k
    reason = sprintf('actuator count m=%d <= active row count k=%d', m, k);
end
end

function u = run_one_allocator(method, y, B, umin, umax, tie_opts)
% 所有算法统一返回限幅后的 u_raw。
[k, m] = size(B);

if nargin < 6 || isempty(tie_opts)
    tie_opts = struct();
end

switch lower(method)
    case 'inv'
        u = pinv(B) * y;

    case 'pca_dir'
        [~, u_tmp, ~, ~] = evalc('DP_LPCA(y, B, umin, umax, 100, tie_opts);');
        u = u_tmp(:);

    case 'pca_dpscaled'
        [~, u_tmp, ~, ~, ~] = evalc('DPscaled_LPCA(y, B, umin, umax, 100, tie_opts);');
        u = u_tmp(:);

    case 'pca_prio'
        m_higher = zeros(k, 1);
        [~, u_tmp, ~, ~] = evalc('DP_LPCA_prio(m_higher, y, B, umin, umax, 100, tie_opts);');
        u = u_tmp(:);

    case 'wls'
        Wv = eye(k);
        Wu = eye(m);
        ud = zeros(m, 1);
        gam = 1e6;
        u0 = (umin + umax) / 2;
        W0 = zeros(m, 1);
        [~, u_tmp, ~, ~] = evalc('wls_alloc(B, y, umin, umax, Wv, Wu, ud, gam, u0, W0, 100);');
        u = u_tmp(:);

    otherwise
        error('未知算法: %s', method);
end

u = clamp_u(u, umin, umax);
end

%% Printing helpers used by main
function print_case_header(case_idx, case_count, B_case, flightData, N)
labels = axis_labels();
fprintf('\n[Case %d/%d] %s\n', case_idx, case_count, B_case.name);
fprintf('  B           : %dx%d\n', size(B_case.B, 1), size(B_case.B, 2));
fprintf('  instances   : %d\n', numel(B_case.inst));
fprintf('  active axes : %s\n', join(labels(B_case.rows), ' '));
fprintf('  unitized    : %s\n', yesno(B_case.unit_allocation_enabled));

if isfield(flightData, 'window')
    fprintf('  window      : %d samples, %.3f to %.3f s\n', ...
        N, flightData.window.t_start_rel_s, flightData.window.t_end_rel_s);
else
    fprintf('  window      : %d samples\n', N);
end

fprintf('\n');
fprintf('  instance  B_raw  B_unit  active axes      u range\n');
fprintf('  --------  -----  ------  ---------------  --------\n');

for i = 1:numel(B_case.inst)
    inst = B_case.inst(i);
    fprintf('  %8d  %dx%-2d  %dx%-3d %-15s  %.3g..%.3g\n', ...
        i - 1, size(inst.B_raw, 1), size(inst.B_raw, 2), ...
        size(inst.B, 1), size(inst.B, 2), ...
        char(join(labels(inst.rows), ' ')), min(inst.umin), max(inst.umax));
end

fprintf('\n');
end

function labels = axis_labels()
labels = ["Mx", "My", "Mz", "Fx", "Fy", "Fz"];
end

function text = yesno(tf)
if tf
    text = 'yes';
else
    text = 'no';
end
end

%% Window and small utilities
function flightData = select_window(flightData, mode, sample_window, time_window)
N = size(flightData.control_sp, 1);
idx = 1:N;

if isfield(flightData, 't') && numel(flightData.t) >= N
    t_rel = flightData.t(:) - flightData.t(1);
else
    t_rel = (0:(N - 1))';
end

switch lower(string(mode))
    case "sample"
        idx = intersect(idx, sample_window);
    case "time"
        idx = intersect(idx, find(t_rel >= time_window(1) & t_rel <= time_window(2)));
    case "all"
    otherwise
        error('未知 WINDOW_MODE=%s。', mode);
end

idx = idx(:)';

if isempty(idx)
    error('窗口没有样本。');
end

fields = {'timestamp_us', 'timestamp_sample_us', 't', 'control_sp', 'u_px4'};

for f = 1:numel(fields)
    name = fields{f};

    if isfield(flightData, name)
        flightData.(name) = crop_rows(flightData.(name), idx);
    end
end

if isfield(flightData, 'setpoint_instances')
    for i = 1:numel(flightData.setpoint_instances)
        for f = 1:numel(fields)
            name = fields{f};

            if isfield(flightData.setpoint_instances(i), name)
                flightData.setpoint_instances(i).(name) = ...
                    crop_rows(flightData.setpoint_instances(i).(name), idx);
            end
        end
    end
end

flightData.window = struct();
flightData.window.full_sample_count = N;
flightData.window.selected_sample_count = numel(idx);
flightData.window.sample_indices = idx;
flightData.window.t_start_rel_s = t_rel(idx(1));
flightData.window.t_end_rel_s = t_rel(idx(end));
end

function value = crop_rows(value, idx)
if (isnumeric(value) || islogical(value)) && ismatrix(value) && size(value, 1) >= max(idx)
    value = value(idx, :);
end
end

function t = time_vector(flightData, N)
if isfield(flightData, 'dt_mean_s') && isfinite(flightData.dt_mean_s)
    dt = flightData.dt_mean_s;
else
    dt = 0.01;
end

t = (0:(N - 1))' * dt;
end

function B = force_6_rows(B)
B = double(B);
B(~isfinite(B)) = 0;

if size(B, 1) < 6
    B(end+1:6, :) = 0;
elseif size(B, 1) > 6
    B = B(1:6, :);
end
end

function Y = force_6_cols(Y)
Y = double(Y);
Y(~isfinite(Y)) = 0;

if size(Y, 2) < 6
    Y(:, end+1:6) = 0;
elseif size(Y, 2) > 6
    Y = Y(:, 1:6);
end
end

function u = clamp_u(u, umin, umax)
u = min(max(u(:), umin(:)), umax(:));
end

function [rmse_u, max_u] = compare_with_px4_output(u, flightData)
% 使用 parser 已经插值到 command sample 的 PX4 actuator 输出。
u_ref = flightData.u_px4;

if isempty(u_ref) || size(u_ref, 2) ~= size(u, 2)
    rmse_u = nan;
    max_u = nan;
    return;
end

N = min(size(u, 1), size(u_ref, 1));
du = u(1:N, :) - u_ref(1:N, :);
v = du(isfinite(du));

if isempty(v)
    rmse_u = nan;
    max_u = nan;
else
    rmse_u = sqrt(mean(v.^2));
    max_u = max(abs(v));
end
end
