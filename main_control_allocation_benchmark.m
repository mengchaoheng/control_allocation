%% Minimal control-allocation benchmark
% 目标很窄：读取当前 PX4 ULog 里的控制分配输入，给一个或多个手写 B case，
% 遍历不同分配算法，比较运行时间、分配误差 y-B*u、以及 u 的差异。
%
% 主流程只有四步：
%   1. 从 PX4 ULog 读取控制分配数据：
%        输入 topic:
%          vehicle_torque_setpoint instance i -> [Mx My Mz]
%          vehicle_thrust_setpoint instance i -> [Fx Fy Fz]
%        parser 会拼成每个 instance 一个 N x 6 的 control_sp：
%          [Mx My Mz Fx Fy Fz]
%        输出参考 topic:
%          actuator_motors.control + actuator_servos.control
%        parser 会插值到输入采样点，保存为 flightData.u_px4。
%   2. 用户在本文件里手写一个或多个 B case：
%        每个 B 实例必须是 6 x n，行顺序 [Mx My Mz Fx Fy Fz]，
%        列顺序就是执行器顺序。
%   3. 预处理 B 和输入：
%        可选 px4_normalize_B() 单位化；
%        找 B 的非零行 rows；
%        每个算法实际只取自己需要的 B/y。
%   4. 对每个 B case、每个算法、每个样本、每个实例分配，最后拼成一个 u。
%
% 结果比较口径：
%   - 主基准固定为“当前 B 上离线重算的 inv”。这样即使日志机型和这里手写 B
%     不同，也仍然有同维度、同 B 的 u 基准。
%   - PX4 日志输出 u_px4 只在列数和当前 B 的执行器数完全一致时比较；
%     如果列数不一致，说明执行器定义不同，只使用日志输入 y 做算法测试。

clc; clear all; close all;

%% 1. 路径和 PX4 日志
tool_dir = fileparts(mfilename('fullpath'));
% 离线测试脚本就是算法库入口，直接加入当前目录下所有已知算法/辅助函数。
% 新增算法文件后通常只需要改 METHODS_TO_RUN。
addpath(genpath(tool_dir));

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

% 样本选择只保留一个简单入口：
%   []      = 全部样本；
%   1:200   = 前 200 个样本；
%   500:900 = 指定样本段。
% 如果要按时间或状态截取，建议在 parser 或日志前处理脚本里处理。
SAMPLE_RANGE = [];

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

case_cfg = struct();
case_cfg.name = "df4_single";
case_cfg.B = [B0 B1];
case_cfg.umin = [0; -ones(4, 1)];
case_cfg.umax = [1;  ones(4, 1)];
case_cfg.normalize_rpy = true;
B_CASES{end+1} = case_cfg;

% SHC09 split 示例：
%   instance 0: motor，主要产生 Fz；
%   instance 1: 六个舵面，产生 Mx/My/Mz。
B0 = [
      0
      0
      0
      0
      0
    -6.5
];

shc09_I = [0.0438 0.0436 0.005006];
shc09_l1 = 0.267;
shc09_l2 = 0.066;
shc09_k_moment = 1.0;
shc09_d = 60*pi/180;
B1_moment = make_ducted_fan_moment_B("hex6", shc09_I, shc09_l1, shc09_l2, shc09_k_moment, shc09_d);
B1 = [B1_moment; zeros(3, size(B1_moment, 2))];

case_cfg = struct();
case_cfg.name = "shc09_split";
case_cfg.B = {B0, B1};
case_cfg.umin = {0, -ones(6, 1)};
case_cfg.umax = {1,  ones(6, 1)};
case_cfg.normalize_rpy = [true true];
B_CASES{end+1} = case_cfg;

% 单实例 B 示例：
case_cfg = struct();
case_cfg.name = "shc09_single";
case_cfg.B = [B0 B1];
case_cfg.umin = [0; -ones(6, 1)];
case_cfg.umax = [1;  ones(6, 1)];
case_cfg.normalize_rpy = true;
B_CASES{end+1} = case_cfg;

% B_CASES 里放什么就跑什么；临时不跑某个 case，直接注释掉对应 block。

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
% 支持的算法名见 run_one_allocator()。这里直接写要跑的算法。
METHODS_TO_RUN = {'inv', 'pca_dir', 'pca_dpscaled', 'wls'};
RESTORING_SELECTION = "restored";    % "raw" / "restored" / "both"
WARMUP_SAMPLE_COUNT = 100;            % MATLAB JIT warmup，不计入算法耗时；设 0 可看冷启动时间。
FALLBACK_METHOD = "inv";
ENABLE_FALLBACK = true;
ENFORCE_U_LIMITS = true;

% 可选：同时调用 alloc_cpp/test/ca_offline_benchmark.cpp 生成 C++ allocator 结果。
% C++ 程序吃同一个 B_case 和同一个日志输入，输出 cpp_pca_dir/cpp_pca_dpscaled。
RUN_CPP_ALLOCATOR = false;
CPP_ALLOCATOR_EXE = fullfile(tool_dir, 'alloc_cpp', 'build', 'ca_offline_benchmark');

% 单位化开关：
%   NORMALIZE_B     : B <- D*B，D 来自 px4_normalize_B.m。
%   NORMALIZE_INPUT : y <- D*y，使用同一个 D。
% 通常：
%   - 日志 y 已经是 PX4 allocator 当时使用的单位：NORMALIZE_INPUT=false。
%   - 物理 B + 物理 y 一起离线测试：NORMALIZE_B=true 且 NORMALIZE_INPUT=true。
NORMALIZE_B = true;
NORMALIZE_INPUT = false;
DEFAULT_NORMALIZE_RPY = true;
ZERO_ROW_TOL = 1e-9;

% main 运行完可以自动打印/画图，不需要手动再运行 print/plot。
RUN_PRINT_AFTER_COMPARE = true;
RUN_PLOT_AFTER_COMPARE = true;
PLOT_SAVE_PNG = false; % 保存会增加运行时间
PLOT_SAVE_FIG = false;
PLOT_SHOW_FIGURES = usejava('desktop');  % 桌面 MATLAB 弹图；batch/headless 只保存。

%% 4. 解析或加载 PX4 日志输入
% 第一次运行时 RESULT_DIR 里没有 command MAT，main 会自动解析日志；
% 以后默认直接加载 MAT，避免把 ulog2csv 的时间混进分配算法测试。
%
% COMMAND_MAT 里只要求有 flightData：
%   flightData.control_sp         默认 instance 的 N x 6 输入；
%   flightData.setpoint_instances 每个 PX4 setpoint instance 的 N x 6 输入；
%   flightData.u_px4              PX4 actuator_motors + actuator_servos 输出参考。
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

if ~isfield(S, 'flightData')
    error('COMMAND_MAT must contain flightData: %s', COMMAND_MAT);
end

flightData = select_samples(S.flightData, SAMPLE_RANGE);

% parser 已经把 vehicle_torque_setpoint + vehicle_thrust_setpoint 拼成 control_sp。
% 这里直接取每个 uORB instance 的 N x 6 输入，列顺序固定为 [Mx My Mz Fx Fy Fz]。
if isfield(flightData, 'setpoint_instances') && ~isempty(flightData.setpoint_instances)
    log_Y = cell(1, numel(flightData.setpoint_instances));

    for i = 1:numel(log_Y)
        log_Y{i} = force_6_cols(flightData.setpoint_instances(i).control_sp);
    end
else
    log_Y = {force_6_cols(flightData.control_sp)};
end

%% 5. 预处理 B，然后进入分配循环
if ~isfolder(RESULT_DIR)
    mkdir(RESULT_DIR);
end

methods_to_run = cellstr(METHODS_TO_RUN);
matlab_methods_to_run = methods_to_run(~cellfun(@is_cpp_allocator_method, methods_to_run));
restore_modes = lower(string(RESTORING_SELECTION));
if restore_modes == "both"
    restore_modes = ["raw", "restored"];
end

RESTORING_AVAILABLE = exist('restoring_cpp', 'file') == 2;
allocator_opts = struct( ...
    'enforce_u_limits', ENFORCE_U_LIMITS, ...
    'enable_fallback', ENABLE_FALLBACK, ...
    'fallback_method', char(FALLBACK_METHOD), ...
    'restoring_available', RESTORING_AVAILABLE);
results_cell = cell(1, numel(B_CASES));
benchmark_tic = tic;

fprintf('\nControl allocation offline benchmark\n');
fprintf('  command mat : %s\n', COMMAND_MAT);
fprintf('  input       : %s\n', make_px4_input_summary(flightData, log_Y));
fprintf('  input data  : %d instance(s), %d sample(s)\n', numel(log_Y), size(log_Y{1}, 1));
fprintf('  methods     : %s\n', strjoin(string(methods_to_run), ', '));
fprintf('  restoring   : %s\n', RESTORING_SELECTION);
fprintf('  warmup      : %d sample(s), not counted in per-method timing\n', WARMUP_SAMPLE_COUNT);

for case_idx = 1:numel(B_CASES)
    B_case = prepare_B_case(B_CASES{case_idx}, NORMALIZE_B, NORMALIZE_INPUT, ...
        DEFAULT_NORMALIZE_RPY, ZERO_ROW_TOL);
    case_run = prepare_case_run(B_case, log_Y, flightData, NORMALIZE_INPUT);
    cpp_process_wall_s = nan;
    alg_cell = cell(1, numel(matlab_methods_to_run) * numel(restore_modes));
    alg_count = 0;

    print_case_header(case_idx, numel(B_CASES), B_case, flightData, case_run.N);

    % 核心测试循环：一个 B case 下面，按算法和 raw/restored 模式逐个跑完整飞行序列。
    for method_idx = 1:numel(matlab_methods_to_run)
        for restore_idx = 1:numel(restore_modes)
            job = struct();
            job.method = lower(char(matlab_methods_to_run{method_idx}));
            job.use_restoring = restore_modes(restore_idx) == "restored";
            job.restoring_mode = char(restore_modes(restore_idx));

            if numel(restore_modes) == 1
                job.name = job.method;
            else
                job.name = sprintf('%s_%s', job.method, job.restoring_mode);
            end

            job_run = prepare_allocator_job(case_run, job, allocator_opts);
            warmup_allocator_job(job_run, WARMUP_SAMPLE_COUNT, allocator_opts);
            alg = run_allocator_series(job_run, case_run, allocator_opts);
            alg_count = alg_count + 1;
            alg_cell{alg_count} = alg;

            fprintf('  %-14s total=%6.4fs avg_alloc=%7.2fus avg_restore=%7.2fus residual_rms=%9.4g fb=%d fail=%d\n', ...
                alg.name, alg.elapsed_s, alg.avg_us_per_sample, alg.avg_restore_us_per_sample, ...
                alg.rmse_residual, alg.fallback_count, alg.fail_count);
        end
    end

    alg_cell = alg_cell(1:alg_count);

    if isempty(alg_cell)
        alg_array = struct([]);
    else
        alg_array = [alg_cell{:}];
    end

    if RUN_CPP_ALLOCATOR || any(cellfun(@is_cpp_allocator_method, methods_to_run))
        [cpp_alg, cpp_process_wall_s] = run_cpp_allocator_case(B_case, case_run.Y_cell, case_run.y_command, flightData, ...
            RESULT_DIR, CPP_ALLOCATOR_EXE, methods_to_run, any(restore_modes == "restored"), case_idx);

        if ~isempty(cpp_alg)
            alg_array = [alg_array cpp_alg];
        end
    end

    [alg_array, u_inv_ref, inv_reference_name] = attach_offline_inv_reference(alg_array);
    results_cell{case_idx} = make_case_result(B_case, case_run, alg_array, ...
        u_inv_ref, inv_reference_name, cpp_process_wall_s);
end

results = [results_cell{:}];
benchmark_elapsed_s = toc(benchmark_tic);

meta = struct();
meta.MODEL_NAME = flightData.model;
meta.COMMAND_MAT = COMMAND_MAT;
meta.LOG_PATH = LOG_PATH;
meta.FORCE_REPARSE = FORCE_REPARSE;
meta.SAMPLE_RANGE = SAMPLE_RANGE;
meta.METHODS_TO_RUN = methods_to_run;
meta.RESTORING_SELECTION = RESTORING_SELECTION;
meta.RESTORING_AVAILABLE = RESTORING_AVAILABLE;
meta.WARMUP_SAMPLE_COUNT = WARMUP_SAMPLE_COUNT;
meta.NORMALIZE_B = NORMALIZE_B;
meta.NORMALIZE_INPUT = NORMALIZE_INPUT;
meta.DEFAULT_NORMALIZE_RPY = DEFAULT_NORMALIZE_RPY;
meta.RUN_CPP_ALLOCATOR = RUN_CPP_ALLOCATOR;
meta.CPP_ALLOCATOR_EXE = CPP_ALLOCATOR_EXE;
meta.benchmark_elapsed_s = benchmark_elapsed_s;
meta.inv_reference_note = 'Primary reference is offline inv recomputed on each current B case.';
meta.reference_note = 'flightData.u_px4 is used only when output column count exactly matches current B.';
meta.px4_reference_note = meta.reference_note;

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

%% B preprocessing
function B_case = prepare_B_case(cfg, normalize_B, normalize_input, default_normalize_rpy, zero_tol)
% 把用户输入的单实例/多实例 B 统一成 6 x n cell，并为每个实例保存：
%   B_raw : 原始 B；
%   B     : 根据 NORMALIZE_B 决定是否单位化后的 B；
%   D     : px4_normalize_B 给出的输入/B 行缩放矩阵；
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
    'D', [], 'scale', [], 'rows', [], 'umin', [], 'umax', [], 'cols', [], 'm', 0), 1, num_inst);

total_cols = 0;

for i = 1:num_inst
    total_cols = total_cols + size(force_6_rows(B_cell{i}), 2);
end

B_all = zeros(6, total_cols);
umin_all = zeros(total_cols, 1);
umax_all = zeros(total_cols, 1);
next_col = 1;

for i = 1:num_inst
    B_raw = force_6_rows(B_cell{i});
    umin = double(umin_cell{i}(:));
    umax = double(umax_cell{i}(:));

    if size(B_raw, 2) ~= numel(umin) || size(B_raw, 2) ~= numel(umax)
        error('%s instance %d: B 列数必须等于 umin/umax 长度。', cfg.name, i - 1);
    end

    normalize_rpy = case_normalize_rpy(cfg, i, default_normalize_rpy);

    if normalize_B || normalize_input
        [D, B_unit, scale] = compute_B_normalization(B_raw, normalize_rpy);
    else
        D = eye(6);
        B_unit = B_raw;
        scale = ones(6, 1);
    end

    if normalize_B
        B = B_unit;
    else
        B = B_raw;
    end

    rows = find(any(abs(B) > zero_tol, 2));
    cols = next_col:(next_col + size(B, 2) - 1);
    next_col = next_col + size(B, 2);

    inst(i).B_raw = B_raw;
    inst(i).B = B;
    inst(i).B_eff = B(rows, :);
    inst(i).D = D;
    inst(i).scale = scale;
    inst(i).rows = rows(:)';
    inst(i).umin = umin;
    inst(i).umax = umax;
    inst(i).cols = cols;
    inst(i).m = size(B, 2);

    B_all(:, cols) = B;
    umin_all(cols) = umin;
    umax_all(cols) = umax;
end

B_case = struct();
B_case.name = char(cfg.name);
B_case.inst = inst;
B_case.B = B_all;
B_case.umin = umin_all;
B_case.umax = umax_all;
B_case.rows = unique([inst.rows]);
B_case.normalize_B_enabled = normalize_B;
B_case.normalize_input_enabled = normalize_input;
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

function [D, B_unit, scale] = compute_B_normalization(B_raw, normalize_rpy)
% 单位化固定使用现有 px4_normalize_B.m，不在这里维护第二套公式。
if exist('px4_normalize_B', 'file') ~= 2
    error('px4_normalize_B.m not found. Check MATLAB path at top of main.');
end

[D, B_unit, ~, scale] = px4_normalize_B(B_raw, normalize_rpy);
B_unit = force_6_rows(B_unit);
end

function text = make_px4_input_summary(flightData, log_Y)
% 打印当前 main 实际使用的数据来源和维数，避免再把输入接口藏起来。
% 数据固定来自 parser：
%   vehicle_torque_setpoint.xyz + vehicle_thrust_setpoint.xyz -> N x 6 control_sp
%   actuator_motors.control + actuator_servos.control -> flightData.u_px4
parts = strings(1, numel(log_Y));

for i = 1:numel(log_Y)
    instance_id = i - 1;

    if isfield(flightData, 'setpoint_instances') && numel(flightData.setpoint_instances) >= i
        instance_id = flightData.setpoint_instances(i).instance;
    end

    parts(i) = sprintf('inst%d:%dx6 [Mx My Mz Fx Fy Fz]', ...
        instance_id, size(log_Y{i}, 1));
end

if isfield(flightData, 'u_px4') && ~isempty(flightData.u_px4)
    ref_text = sprintf('u_px4:%dx%d', size(flightData.u_px4, 1), size(flightData.u_px4, 2));
else
    ref_text = 'u_px4:empty';
end

text = sprintf('PX4 topics vehicle_torque_setpoint + vehicle_thrust_setpoint, %s, %s', ...
    char(strjoin(parts, '; ')), ref_text);
end

function case_run = prepare_case_run(B_case, log_Y, flightData, normalize_input)
% 把日志输入映射到当前 B case：
%   日志实例数 == B 实例数：一一对应；
%   日志多实例、B 单实例：多个 y 相加给一个大 B；
%   日志单实例、B 多实例：force/Fz 给第 0 个 B，torque/Mxyz 给第 1 个 B。
num_log = numel(log_Y);
num_B = numel(B_case.inst);

if num_log == num_B
    Y_cell = log_Y;
elseif num_log > 1 && num_B == 1
    N_sum = min(cellfun(@(Y) size(Y, 1), log_Y));
    Y_sum = zeros(N_sum, 6);

    for i = 1:num_log
        Y_sum = Y_sum + log_Y{i}(1:N_sum, :);
    end

    Y_cell = {Y_sum};
elseif num_log == 1 && num_B > 1
    Y = log_Y{1};
    Y_cell = cell(1, num_B);

    for i = 1:num_B
        Y_cell{i} = zeros(size(Y));
    end

    Y_cell{1}(:, 4:6) = Y(:, 4:6);
    Y_cell{2}(:, 1:3) = Y(:, 1:3);
else
    error('不支持日志实例数 %d -> B 实例数 %d 的输入转换。', num_log, num_B);
end

if normalize_input
    for i = 1:numel(Y_cell)
        Y_cell{i} = Y_cell{i} * B_case.inst(i).D';
    end
end

N = min(cellfun(@(Y) size(Y, 1), Y_cell));
case_run = B_case;
y_command = zeros(N, 6);

for inst_idx = 1:numel(case_run.inst)
    Y_full = Y_cell{inst_idx}(1:N, :);
    case_run.inst(inst_idx).Y_full = Y_full;
    case_run.inst(inst_idx).Y_eff = Y_full(:, case_run.inst(inst_idx).rows);
    Y_cell{inst_idx} = Y_full;
    y_command = y_command + Y_full;
end

case_run.N = N;
case_run.Y_cell = Y_cell;
case_run.y_command = y_command;

if isfield(flightData, 't') && numel(flightData.t) >= N
    case_run.t = flightData.t(1:N);
else
    case_run.t = (0:(N - 1))' * flightData.dt_mean_s;
end

case_run.flightData = flightData;
end

%% Allocation
function job_run = prepare_allocator_job(case_run, job, opts)
% 每个算法对每个实例只准备一次：
%   - PCA/LPCA 类只看有效行 B_eff/Y_eff；
%   - 其他算法直接看完整 6 行 B/Y；
%   - 退化 B 是否需要 fallback 也在这里判断一次。
check_method_available(job.method);

job_run = job;
job_run.N = case_run.N;
job_run.total_u_dim = numel(case_run.umin);
inst_template = struct('B_full', [], 'Y_full', [], 'B_restore', [], ...
    'rows', [], 'cols', [], 'm', 0, ...
    'umin_check', [], 'umax_check', [], ...
    'B_run', [], 'Y_run', [], 'umin', [], 'umax', [], ...
    'fallback_reason', "");
job_run.inst = repmat(inst_template, 1, numel(case_run.inst));

for inst_idx = 1:numel(case_run.inst)
    inst = case_run.inst(inst_idx);
    jr = struct();
    jr.B_full = inst.B;
    jr.Y_full = inst.Y_full;
    jr.B_restore = inst.B_eff;
    jr.rows = inst.rows;
    jr.cols = inst.cols;
    jr.m = inst.m;
    jr.umin_check = inst.umin;
    jr.umax_check = inst.umax;

    if method_uses_effective_rows(job.method)
        jr.B_run = inst.B_eff;
        jr.Y_run = inst.Y_eff;
    else
        jr.B_run = inst.B;
        jr.Y_run = inst.Y_full;
    end

    if opts.enforce_u_limits
        jr.umin = inst.umin;
        jr.umax = inst.umax;
    else
        jr.umin = -inf(size(inst.umin));
        jr.umax = inf(size(inst.umax));
    end

    jr.fallback_reason = direct_fallback_reason(job.method, jr.B_run);
    job_run.inst(inst_idx) = jr;
end
end

function warmup_allocator_job(job_run, warmup_count, opts)
% MATLAB 第一次调用某些 allocator 会有 JIT/函数加载开销。
% warmup 只预热，不保存，不计入 allocator 平均耗时。
if warmup_count <= 0
    return;
end

N = min(warmup_count, job_run.N);

for sample_idx = 1:N
    for inst_idx = 1:numel(job_run.inst)
        jr = job_run.inst(inst_idx);

        if isempty(jr.rows)
            continue;
        end

        try
            y_run = jr.Y_run(sample_idx, :)';
            y_full = jr.Y_full(sample_idx, :)';
            u_raw = run_allocator_sample(job_run.method, y_run, jr.B_run, ...
                y_full, jr.B_full, jr.umin, jr.umax, jr.fallback_reason, opts);
            apply_optional_restoring(jr.B_restore, u_raw, jr.umin, jr.umax, ...
                job_run.use_restoring, opts.restoring_available);
        catch
            % 正式循环会记录失败样本；warmup 不干扰 benchmark 主流程。
        end
    end
end
end

function alg = run_allocator_series(job_run, case_run, opts)
% 真正的分配循环只保留三件事：
%   1. 对每个样本、每个实例求 u；
%   2. 把多实例 u 按执行器列拼起来；
%   3. 计算 y-B*u、运行时间和失败/fallback 计数。
N = case_run.N;
u_all = nan(N, job_run.total_u_dim);
y_achieved = nan(N, 6);
residual = nan(N, 6);

primary_time = 0;
fallback_time = 0;
restore_time = 0;
fail_count = 0;
fallback_count = 0;
alg_tic = tic;

for sample_idx = 1:N
    y_out = zeros(6, 1);
    sample_failed = false;

    for inst_idx = 1:numel(job_run.inst)
        jr = job_run.inst(inst_idx);

        if isempty(jr.rows)
            u = zeros(jr.m, 1);
            did_fallback = false;
        else
            y_run = jr.Y_run(sample_idx, :)';
            y_full = jr.Y_full(sample_idx, :)';

            try
                [u_raw, did_fallback, primary_dt, fallback_dt] = run_allocator_sample( ...
                    job_run.method, y_run, jr.B_run, y_full, jr.B_full, ...
                    jr.umin, jr.umax, jr.fallback_reason, opts);
                primary_time = primary_time + primary_dt;
                fallback_time = fallback_time + fallback_dt;

                [u, restore_dt] = apply_optional_restoring(jr.B_restore, u_raw, ...
                    jr.umin, jr.umax, job_run.use_restoring, opts.restoring_available);
                restore_time = restore_time + restore_dt;
            catch
                u = nan(jr.m, 1);
                did_fallback = false;
                sample_failed = true;
            end
        end

        u_all(sample_idx, jr.cols) = u(:)';
        y_out = y_out + jr.B_full * u;
        fallback_count = fallback_count + double(did_fallback);

        if ~all(isfinite(u)) || any(u < jr.umin_check - 1e-7) || any(u > jr.umax_check + 1e-7)
            sample_failed = true;
        end
    end

    y_achieved(sample_idx, :) = y_out';
    residual(sample_idx, :) = case_run.y_command(sample_idx, :) - y_out';
    fail_count = fail_count + double(sample_failed);
end

allocator_time = primary_time + fallback_time;
elapsed_s = toc(alg_tic);
[residual_rms, residual_max] = finite_rms_max(residual);
[rmse_vs_px4, max_abs_vs_px4] = compare_with_px4_output(u_all, case_run.flightData);

alg = struct( ...
    'name', job_run.name, ...
    'method', job_run.method, ...
    'restoring_mode', job_run.restoring_mode, ...
    'u', u_all, ...
    'y_achieved', y_achieved, ...
    'residual', residual, ...
    'elapsed_s', elapsed_s, ...
    'allocator_s', allocator_time, ...
    'restore_s', restore_time, ...
    'avg_us_per_sample', 1e6 * allocator_time / max(N, 1), ...
    'avg_restore_us_per_sample', 1e6 * restore_time / max(N, 1), ...
    'rmse_residual', residual_rms, ...
    'max_abs_residual', residual_max, ...
    'rmse_vs_inv', nan, ...
    'max_abs_vs_inv', nan, ...
    'rmse_vs_px4', rmse_vs_px4, ...
    'max_abs_vs_px4', max_abs_vs_px4, ...
    'fail_count', fail_count, ...
    'fallback_count', fallback_count);
end

function [u, did_fallback, primary_dt, fallback_dt] = run_allocator_sample( ...
    method, y_run, B_run, y_full, B_full, umin, umax, fallback_reason, opts)
% 单个样本的 allocator 调用。直接 fallback 的条件已经提前算好，
% 这里不再做 rank 等重复判断。
primary_dt = 0;
fallback_dt = 0;
did_fallback = false;

if strlength(fallback_reason) > 0
    if ~opts.enable_fallback
        error('%s cannot run on this B: %s', method, fallback_reason);
    end

    fallback_tic = tic;
    u = run_one_allocator(opts.fallback_method, y_full, B_full, umin, umax);
    fallback_dt = toc(fallback_tic);
    did_fallback = true;
    return;
end

try
    primary_tic = tic;
    u = run_one_allocator(method, y_run, B_run, umin, umax);
    primary_dt = toc(primary_tic);
catch ME
    primary_dt = toc(primary_tic);

    if ~opts.enable_fallback
        rethrow(ME);
    end

    fallback_tic = tic;
    u = run_one_allocator(opts.fallback_method, y_full, B_full, umin, umax);
    fallback_dt = toc(fallback_tic);
    did_fallback = true;
end
end

function [u, elapsed_s] = apply_optional_restoring(B, u_raw, umin, umax, use_restoring, restoring_available)
% raw/restored 只在这里切换，方便比较恢复算法前后的结果。
if use_restoring && restoring_available && ~isempty(B)
    restore_tic = tic;
    u = restoring_cpp(B, u_raw(:), umin, umax);
    elapsed_s = toc(restore_tic);
else
    u = u_raw(:);
    elapsed_s = 0;
end
end

function reason = direct_fallback_reason(method, B)
% PCA/LPCA 需要有效行空间里有自由度。k<=1 或 m<=k 时，
% 这个实例没有可优化零空间，直接用 inv 作为 fallback。
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

function tf = method_uses_effective_rows(method)
% 目前只有 PCA/LPCA 这类算法吃有效行 B_eff；其他算法保留 6 x n。
tf = any(strcmpi(method, {'pca_dir', 'pca_dpscaled', 'pca_prio'}));
end

function tf = is_cpp_allocator_method(method)
tf = startsWith(lower(string(method)), "cpp_");
end

function check_method_available(method)
% 缺算法文件时直接报错，不在样本循环里反复判断。
switch lower(method)
    case 'pca_dir'
        required = {'DP_LPCA'};
    case 'pca_dpscaled'
        required = {'DPscaled_LPCA'};
    case 'pca_prio'
        required = {'DP_LPCA_prio'};
    case {'lib_lpwrap_db', 'lib_lpwrap_dbinf', 'lib_lpwrap_dir', 'lib_lpwrap_dpscaled', ...
            'lib_lpwrap_mo', 'lib_lpwrap_sb', 'lib_lpwrap_incre'}
        required = {'LPwrap'};
    case {'lib_lpwrap_par_db', 'lib_lpwrap_par_dbinf', 'lib_lpwrap_par_dir', ...
            'lib_lpwrap_par_dpscaled', 'lib_lpwrap_par_mo', 'lib_lpwrap_par_sb'}
        required = {'LPwrap_par'};
    case 'lib_cgiwrap'
        required = {'CGIwrap'};
    case 'lib_dawrap'
        required = {'DAwrap'};
    case 'lib_vjawrap'
        required = {'VJAwrap'};
    case 'wls'
        required = {'wls_alloc'};
    case 'wls_gen'
        required = {'wls_alloc_gen'};
    case 'dir_linprog'
        required = {'dir_alloc_linprog'};
    case 'dir_linprog_re'
        required = {'dir_alloc_linprog_re'};
    case 'dir_linprog_re_bound'
        required = {'dir_alloc_linprog_re_bound'};
    case 'use_lp_lib'
        required = {'use_LP_lib'};
    case 'allocator_dir_lpwrap_4'
        required = {'allocator_dir_LPwrap_4'};
    otherwise
        required = {};
end

for i = 1:numel(required)
    if exist(required{i}, 'file') ~= 2
        error('%s.m not found for method %s.', required{i}, method);
    end
end
end

function u = run_one_allocator(method, y, B, umin, umax)
% 算法分发只接受一个样本的 y/B/限幅，并统一返回限幅后的 u_raw。
[k, m] = size(B);

switch lower(method)
    case 'inv'
        u = pinv(B) * y;

    case 'pca_dir'
        [~, u_tmp, ~, ~] = evalc('DP_LPCA(y, B, umin, umax, 100);');
        u = u_tmp(:);

    case 'pca_dpscaled'
        [~, u_tmp, ~, ~, ~] = evalc('DPscaled_LPCA(y, B, umin, umax, 100);');
        u = u_tmp(:);

    case 'pca_prio'
        m_higher = zeros(k, 1);
        [~, u_tmp, ~, ~] = evalc('DP_LPCA_prio(m_higher, y, B, umin, umax, 100);');
        u = u_tmp(:);

    case 'lib_lpwrap_db'
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 0));');
        u = u_tmp(:);

    case 'lib_lpwrap_dbinf'
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 1));');
        u = u_tmp(:);

    case 'lib_lpwrap_dir'
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 2));');
        u = u_tmp(:);

    case 'lib_lpwrap_dpscaled'
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 3));');
        u = u_tmp(:);

    case 'lib_lpwrap_mo'
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 4));');
        u = u_tmp(:);

    case 'lib_lpwrap_sb'
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 5));');
        u = u_tmp(:);

    case 'lib_lpwrap_par_db'
        [~, u_tmp] = evalc('LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 0), y, m);');
        u = u_tmp(:);

    case 'lib_lpwrap_par_dbinf'
        [~, u_tmp] = evalc('LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 1), y, m);');
        u = u_tmp(:);

    case 'lib_lpwrap_par_dir'
        [~, u_tmp] = evalc('LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 2), y, m);');
        u = u_tmp(:);

    case 'lib_lpwrap_par_dpscaled'
        [~, u_tmp] = evalc('LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 3), y, m);');
        u = u_tmp(:);

    case 'lib_lpwrap_par_mo'
        [~, u_tmp] = evalc('LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 4), y, m);');
        u = u_tmp(:);

    case 'lib_lpwrap_par_sb'
        [~, u_tmp] = evalc('LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 5), y, m);');
        u = u_tmp(:);

    case 'lib_lpwrap_incre'
        u0 = zeros(m, 1);
        [~, u_tmp] = evalc('LPwrap(make_lpwrap_in_mat(B, y, umin - u0, umax - u0, 2));');
        u = u_tmp(:) + u0;

    case 'lib_cgiwrap'
        [~, u_tmp] = evalc('CGIwrap(make_book_wrapper_in_mat(B, y, umin, umax));');
        u = u_tmp(:);

    case 'lib_dawrap'
        [~, u_tmp] = evalc('DAwrap(make_book_wrapper_in_mat(B, y, umin, umax));');
        u = u_tmp(:);

    case 'lib_vjawrap'
        [~, u_tmp] = evalc('VJAwrap(make_book_wrapper_in_mat(B, y, umin, umax));');
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

    case 'wls_gen'
        Wv = eye(k);
        Wu = eye(m);
        ud = zeros(m, 1);
        gam = 1e6;
        u0 = zeros(m, 1);
        W0 = zeros(m, 1);
        [~, u_tmp] = evalc('wls_alloc_gen(B, y, umin, umax, Wv, Wu, ud, gam, u0, W0, 100, m);');
        u = u_tmp(:);

    case 'dir_linprog'
        [~, u_tmp, ~] = evalc('dir_alloc_linprog(B, y, umin, umax, 1e4);');
        u = u_tmp(:);

    case 'dir_linprog_re'
        [~, u_tmp, ~] = evalc('dir_alloc_linprog_re(B, y, umin, umax);');
        u = u_tmp(:);

    case 'dir_linprog_re_bound'
        [~, u_tmp, ~] = evalc('dir_alloc_linprog_re_bound(B, y, umin, umax, 1e4);');
        u = u_tmp(:);

    case 'use_lp_lib'
        [~, u_tmp, ~] = evalc('use_LP_lib(B, y, umin, umax);');
        u = u_tmp(:);

    case 'allocator_dir_lpwrap_4'
        if m ~= 4
            error('allocator_dir_LPwrap_4 only supports m=4, current m=%d.', m);
        end
        [~, u_tmp, ~, ~] = evalc('allocator_dir_LPwrap_4(single(B), single(y), single(umin), single(umax));');
        u = double(u_tmp(:));

    otherwise
        error('未知算法: %s', method);
end

u = clamp_u(u, umin, umax);
end

function IN_MAT = make_lpwrap_in_mat(B, y, umin, umax, lp_method)
% control_allocation_lib/LPwrap 的固定输入格式：[B y; umin 0; umax 0; active lp_method]。
[~, m] = size(B);
active_effectors = ones(1, m);
IN_MAT = [B                y(:);
          umin(:)'         0;
          umax(:)'         0;
          active_effectors lp_method];
end

function IN_MAT = make_book_wrapper_in_mat(B, y, umin, umax)
% CGIwrap/DAwrap/VJAwrap 使用同样矩阵格式，但最后一个标量不是 LPmethod。
IN_MAT = make_lpwrap_in_mat(B, y, umin, umax, 0);
end

function [rms_value, max_value] = finite_rms_max(X)
v = X(isfinite(X));

if isempty(v)
    rms_value = nan;
    max_value = nan;
else
    rms_value = sqrt(mean(v.^2));
    max_value = max(abs(v));
end
end

function result = make_case_result(B_case, case_run, alg_array, u_inv_ref, inv_reference_name, cpp_process_wall_s)
% 保存结果时只保留后续打印/绘图/二次分析真正使用的内容。
result = struct( ...
    'name', B_case.name, ...
    'B', B_case.B, ...
    'umin', B_case.umin, ...
    'umax', B_case.umax, ...
    'rows', B_case.rows, ...
    'cpp_process_wall_s', cpp_process_wall_s, ...
    't', case_run.t, ...
    'control_sp', case_run.y_command, ...
    'u_inv_ref', u_inv_ref, ...
    'inv_reference_name', inv_reference_name, ...
    'alg', alg_array);
end

function cpp_methods = cpp_methods_from_selection(methods_to_run)
% C++ 桥接程序输出文件仍用 pca_dir/pca_dpscaled 命名；
% 用户选择 cpp_dir/cpp_dpscaled 时在这里映射过去。
cpp_methods = strings(numel(methods_to_run), 1);
count = 0;

for i = 1:numel(methods_to_run)
    method = lower(string(methods_to_run{i}));

    switch method
        case {"pca_dir", "cpp_dir"}
            count = count + 1;
            cpp_methods(count) = "pca_dir";
        case {"pca_dpscaled", "cpp_dpscaled"}
            count = count + 1;
            cpp_methods(count) = "pca_dpscaled";
    end
end

cpp_methods = cpp_methods(1:count);
cpp_methods = unique(cpp_methods, 'stable');
end

%% C++ offline allocator bridge
function [cpp_alg, cpp_process_wall_s] = run_cpp_allocator_case(B_case, Y_cell, y_command, flightData, ...
    result_dir, cpp_exe, methods_to_run, cpp_use_restoring, case_idx)
% 把当前 MATLAB main 中已经准备好的 B/Y 交给 alloc_cpp/test/ca_offline_benchmark.cpp。
% C++ 输出会被读回成和 MATLAB allocator 一样的 alg 结构，因此 print/plot 不需要特殊分支。
cpp_alg = struct([]);
cpp_process_wall_s = nan;
cpp_methods = cpp_methods_from_selection(methods_to_run);

if isempty(cpp_methods)
    return;
end

if ~isfile(cpp_exe)
    warning(['C++ allocator executable not found:\n  %s\n', ...
        'Build it once in alloc_cpp/build, then rerun main.'], cpp_exe);
    return;
end

N = size(y_command, 1);
case_stem = sprintf('%02d_%s', case_idx, safe_file_name(B_case.name));
cpp_input_dir = fullfile(result_dir, 'cpp_inputs', case_stem);
cpp_output_dir = fullfile(result_dir, 'cpp_outputs', case_stem);

if ~isfolder(cpp_input_dir)
    mkdir(cpp_input_dir);
end

if ~isfolder(cpp_output_dir)
    mkdir(cpp_output_dir);
end

write_cpp_case_inputs(cpp_input_dir, B_case, Y_cell, N, cpp_use_restoring);

cmd = sprintf('"%s" "%s" "%s"', cpp_exe, cpp_input_dir, cpp_output_dir);
cpp_process_tic = tic;
[status, cmdout] = system(cmd);
cpp_process_wall_s = toc(cpp_process_tic);

if status ~= 0
    warning('C++ allocator failed:\n%s', cmdout);
    return;
end

fprintf('\nC++ allocator bridge:\n');
fprintf('  executable : %s\n', cpp_exe);
fprintf('  input dir  : %s\n', cpp_input_dir);
fprintf('  output dir : %s\n', cpp_output_dir);
fprintf('  process wall: %.4fs (includes process start + CSV I/O)\n', cpp_process_wall_s);

cpp_alg_cell = cell(1, numel(cpp_methods));

for i = 1:numel(cpp_methods)
    method = char(cpp_methods(i));
    cpp_alg_cell{i} = read_cpp_algorithm_result(cpp_output_dir, method, ...
        y_command, flightData, N, cpp_use_restoring);

    fprintf('  %-14s total=%6.4fs avg_alloc=%7.2fus residual_rms=%9.4g fb=%d fail=%d\n', ...
        cpp_alg_cell{i}.name, cpp_alg_cell{i}.elapsed_s, ...
        cpp_alg_cell{i}.avg_us_per_sample, cpp_alg_cell{i}.rmse_residual, ...
        cpp_alg_cell{i}.fallback_count, cpp_alg_cell{i}.fail_count);
end

cpp_alg = [cpp_alg_cell{:}];
end

function write_cpp_case_inputs(input_dir, B_case, Y_cell, N, use_restoring)
% C++ 侧只需要每个 allocation instance 的 B、输入 Y、限幅。
% 文件名保持机械化，便于别的测试程序复用：
%   inst_0_B.csv / inst_0_Y.csv / inst_0_umin.csv / inst_0_umax.csv
case_meta = [N, numel(B_case.umin), numel(B_case.inst), double(use_restoring)];
writematrix(case_meta, fullfile(input_dir, 'case_meta.csv'));

for inst_idx = 1:numel(B_case.inst)
    inst = B_case.inst(inst_idx);
    prefix = fullfile(input_dir, sprintf('inst_%d', inst_idx - 1));
    writematrix(inst.B, [prefix '_B.csv']);
    writematrix(Y_cell{inst_idx}(1:N, :), [prefix '_Y.csv']);
    writematrix(inst.umin(:)', [prefix '_umin.csv']);
    writematrix(inst.umax(:)', [prefix '_umax.csv']);
end
end

function alg = read_cpp_algorithm_result(output_dir, method, y_command, flightData, N, use_restoring)
% C++ 输出 CSV 后，MATLAB 只负责读回来并计算同样的统计量。
u = readmatrix(fullfile(output_dir, [method '_u.csv']));
y_achieved = readmatrix(fullfile(output_dir, [method '_y_achieved.csv']));
residual = readmatrix(fullfile(output_dir, [method '_residual.csv']));
timing = readmatrix(fullfile(output_dir, [method '_timing.csv']));

u = fixed_rows(u, N);
y_achieved = fixed_rows(y_achieved, N);
residual = fixed_rows(residual, N);
residual_values = residual(isfinite(residual));
[rmse_vs_px4, max_abs_vs_px4] = compare_with_px4_output(u, flightData);

if use_restoring
    restoring_mode = 'restored';
else
    restoring_mode = 'raw';
end

alg = struct( ...
    'name', ['cpp_' method], ...
    'method', ['cpp_' method], ...
    'restoring_mode', restoring_mode, ...
    'u', u, ...
    'y_achieved', y_achieved, ...
    'residual', residual, ...
    'elapsed_s', timing(1), ...
    'allocator_s', timing(2), ...
    'restore_s', timing(3), ...
    'avg_us_per_sample', 1e6 * timing(2) / max(N, 1), ...
    'avg_restore_us_per_sample', 1e6 * timing(3) / max(N, 1), ...
    'rmse_residual', sqrt(mean(residual_values.^2)), ...
    'max_abs_residual', max(abs(residual_values)), ...
    'rmse_vs_inv', nan, ...
    'max_abs_vs_inv', nan, ...
    'rmse_vs_px4', rmse_vs_px4, ...
    'max_abs_vs_px4', max_abs_vs_px4, ...
    'fail_count', round(timing(4)), ...
    'fallback_count', round(timing(5)));
end

function [alg_array, u_inv_ref, inv_reference_name] = attach_offline_inv_reference(alg_array)
% 对每个 B case，主参考都用同一个 B 上离线重算出来的 inv。
% 这样日志机型、日志 actuator 维数、当前手写 B 不一致时，仍然可以比较不同算法的 u。
u_inv_ref = [];
inv_reference_name = "";

if isempty(alg_array)
    return;
end

method_names = string({alg_array.method});
restoring_modes = string({alg_array.restoring_mode});
inv_candidates = find(method_names == "inv");

if isempty(inv_candidates)
    warning('No inv result found. u-u_inv reference metrics are disabled for this B case.');
    return;
end

% 如果同时跑 raw/restored，优先用 restored inv，因为多数情况下最终分配结果是 restoring 后的 u。
preferred = inv_candidates(restoring_modes(inv_candidates) == "restored");

if isempty(preferred)
    preferred = inv_candidates(1);
else
    preferred = preferred(1);
end

u_inv_ref = alg_array(preferred).u;
inv_reference_name = string(alg_array(preferred).name);

for i = 1:numel(alg_array)
    [alg_array(i).rmse_vs_inv, alg_array(i).max_abs_vs_inv] = ...
        compare_u_reference(alg_array(i).u, u_inv_ref);
end
end

function X = fixed_rows(X, N)
% readmatrix 读单列/单行 CSV 时形状可能退化；这里只保证样本行数一致。
if size(X, 1) < N
    X(end+1:N, :) = nan;
elseif size(X, 1) > N
    X = X(1:N, :);
end
end

function name = safe_file_name(name)
% 和 plot 文件保持同一个文件名规则。
name = regexprep(char(name), '[^a-zA-Z0-9_\-]+', '_');
name = regexprep(name, '_+', '_');
name = strip(string(name), '_');

if strlength(name) == 0
    name = "case";
end

name = char(name);
end

%% Printing helpers used by main
function print_case_header(case_idx, case_count, B_case, flightData, N)
labels = axis_labels();
fprintf('\n[Case %d/%d] %s\n', case_idx, case_count, B_case.name);
fprintf('  B           : %dx%d\n', size(B_case.B, 1), size(B_case.B, 2));
fprintf('  instances   : %d\n', numel(B_case.inst));
fprintf('  active axes : %s\n', join(labels(B_case.rows), ' '));
fprintf('  normalize B : %s\n', yesno(B_case.normalize_B_enabled));
fprintf('  normalize y : %s\n', yesno(B_case.normalize_input_enabled));

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

%% Sample selection and small utilities
function flightData = select_samples(flightData, sample_range)
% main 只按样本编号截取。复杂的数据选择，例如按飞行状态/时间段截取，
% 建议在 parser 或日志预处理脚本里完成。
N = size(flightData.control_sp, 1);
idx = 1:N;

if isfield(flightData, 't') && numel(flightData.t) >= N
    t_rel = flightData.t(:) - flightData.t(1);
else
    t_rel = (0:(N - 1))';
end

if nargin >= 2 && ~isempty(sample_range)
    idx = intersect(idx, sample_range);
end

idx = idx(:)';

if isempty(idx)
    error('SAMPLE_RANGE 没有选中任何样本。');
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

function B = force_6_rows(B)
B = double(B);
B(~isfinite(B)) = 0;

if size(B, 1) < 6
    B(end+1:6, :) = 0;
elseif size(B, 1) > 6
    error('B must have at most 6 rows in [Mx My Mz Fx Fy Fz]. Current B is %dx%d.', size(B, 1), size(B, 2));
end
end

function Y = force_6_cols(Y)
Y = double(Y);
Y(~isfinite(Y)) = 0;

if size(Y, 2) < 6
    Y(:, end+1:6) = 0;
elseif size(Y, 2) > 6
    error('Control input must have at most 6 columns in [Mx My Mz Fx Fy Fz]. Current Y is %dx%d.', size(Y, 1), size(Y, 2));
end
end

function u = clamp_u(u, umin, umax)
u = min(max(u(:), umin(:)), umax(:));
end

function [rmse_u, max_u] = compare_with_px4_output(u, flightData)
% 使用 parser 已经插值到 command sample 的 PX4 actuator 输出。
u_ref = flightData.u_px4;

% 这里必须要求列数完全相等。只截断/补齐会把“不同机型/不同 actuator 定义”的
% 日志输出误当成同一个分配问题的参考。
if isempty(u_ref) || size(u_ref, 2) ~= size(u, 2)
    rmse_u = nan;
    max_u = nan;
    return;
end

[rmse_u, max_u] = compare_u_reference(u, u_ref);
end

function [rmse_u, max_u] = compare_u_reference(u, u_ref)
% 同维度 u 参考比较工具。inv 参考和 PX4 日志参考都走这个函数。
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
