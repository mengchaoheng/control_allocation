clear all;
close all;
clc;
addpath(genpath(pwd))

%% setup aircraft and load input data
% Original derivation notes kept here intentionally. They are useful when
% switching B definitions or comparing direct inverse style allocators.
%
% 行满秩时pinv(B)=The Moore–Penrose Pseudo-inverse
    % // If B is full raw rank, The Moore–Penrose Pseudo-inverse B^+= B^T (B B^T)^{-1},since
    % // B=K*P, K=I\diag([l1 l1 l2])k, P=[-1     0     1     0; 0    -1     0     1; 1     1     1     1];   K=diag([ k_omega2force*l1/I_x  k_omega2force*l1/I_y  k_omega2force*l2/I_z  ])
    %     // B^{\dagger} = P^\top K K^{-1} (P P^\top)^{-1} K^{-1} = P^\top (P P^\top)^{-1} K^{-1}=P^{\dagger} K^{-1}
    % // P^{\dagger}=[-0.5000   -0.0000    0.2500;0   -0.5000    0.2500;0.5000   -0.0000    0.2500;0    0.5000    0.2500]
    % // B^{\dagger}=P^{\dagger} K^{-1}=[-0.5000   -0.0000    0.2500;0   -0.5000    0.2500;0.5000   -0.0000    0.2500;0    0.5000    0.2500]*diag([ I_x/(k_omega2force*l1)  I_y/(k_omega2force*l1)  I_z/(k_omega2force*l2)  ])

load 'input.mat'; % provides v, len_command_px4, controls_delta_t_s, ...
if isfile('input.csv')
    v = readmatrix('input.csv')';
end

[~, N] = size(v);
dt = mean(controls_delta_t_s);
t_full = 0:dt:dt*(len_command_px4-1);

%% Test configuration
test_time_window_s = [0 4];  % 测试时间窗口；[] = 全程，示例：[20 30] 只跑 20-30 秒。[24.8 25.4]
use_restoring = 'both';            % true=只看 restoring；false=只看 raw allocation；'both'=raw/restoring 各画一张。

allocation_method_selection = {'pca_dir','cpp_dir'};  % 参与测试的分配算法；'all' = 全部，示例：[1 2 5] 或 {'inv','pca_dir','wls_gen','cpp_dir'}。
% 常用名含义：
%   inv: pinv(B)*y 后限幅，再按 use_restoring 决定是否 restoring_cpp。
%   pca_dir: /PCA/DP_LPCA 后按 use_restoring 决定是否 restoring_cpp。
%   pca_dpscaled: /PCA/DPscaled_LPCA 后按 use_restoring 决定是否 restoring_cpp。
%   pca_prio: /PCA/DP_LPCA_prio 后按 use_restoring 决定是否 restoring_cpp。
%   cpp_dir/cpp_dpscaled/cpp_prio: 分别对应 C++ DP_LPCA / DPscaled_LPCA / DP_LPCA_prio 输出。
%   lib_lpwrap_* / lib_lpwrap_par_*: 分配库 aircraft-control-allocation-book-simulation 的 LPwrap.m / LPwrap_par.m，使用 LPmethod=0..5。
%       LPmethod=2 使用 lib_lpwrap_dir / lib_lpwrap_par_dir。
%       LPmethod=3 使用 lib_lpwrap_dpscaled / lib_lpwrap_par_dpscaled。
%   lib_cgiwrap / lib_dawrap / lib_vjawrap: 分配库里的独立包装算法，输入矩阵最后一格按原接口填 0。
%   wls / wls_gen: QCAT WLS 和生成版 WLS。
%   dir_linprog*: linprog 形式的直接分配。
%   use_lp_lib: reformula_LP/use_LP_lib。
%   allocator_dir_lpwrap_4: 只支持 4 执行器的 C 风格 DPscaled 包装，用在非 4 维 B 会报告 fail。

aircraft_case_selection = {'15006_SHC09'};  % 参与测试的 B/case；'all' = 全部，示例：[2] 或 {'15006_SHC09'}。
% case 编号/正式名：
%   1: '15003_ductedfan4'
%   2: '15006_SHC09'
%   3: '15008_SHW09_vtol_full8'
%   4: '15008_SHW09_vtol_reduced6'

reports_to_run = {'allocator_summary','method_diff'};  % 输出报告；'all' = 全部，示例：{'allocator_summary','case_diff'}。
%   allocator_summary: 每个 B 下，每个算法的 ok/fail、误差、饱和情况。
%   method_diff: 相同 B，不同分配算法之间对比，并画同 B 不同算法的 u 曲线叠图。
%   case_diff: 相同分配算法，不同 B/case 之间对比，并画同算法不同 B 的 u 曲线叠图。
%   reachable_set: 使用 QCAT/vview 画可达集。

allocator_sample_print_count = 0;  % 逐样本打印数量；0 = 不打印，示例：3 表示每个算法打印前 3 个样本的 y/u/Bu。

case_comparison_pairs = {  % case_diff 对比组合；每行 {case A, case B}，建议使用上面的正式名。
    '15008_SHW09_vtol_full8', '15008_SHW09_vtol_reduced6'
};

allocation_method_catalog = get_allocation_method_catalog();
allocation_methods_to_run = resolve_allocation_method_selection(allocation_method_selection, allocation_method_catalog);

if isempty(test_time_window_s)
    test_sample_indices = 1:len_command_px4;
else
    test_sample_indices = find(t_full >= test_time_window_s(1) & t_full <= test_time_window_s(2));
end
if isempty(test_sample_indices)
    error('test:EmptyTimeWindow', 'No samples found in test_time_window_s = [%g %g].', test_time_window_s(1), test_time_window_s(2));
end

v_test = v(:, test_sample_indices);
command_px4 = v_test;
len_command_px4 = size(v_test, 2);
t = t_full(test_sample_indices);

%% Aircraft models
aircraft_cases = {make_aircraft_4(), make_aircraft_6()};
aircraft_cases{end+1} = make_shw09_vtol_mc_full8();
aircraft_cases{end+1} = make_shw09_vtol_mc_reduced6();
% To test another B, append a case here:
% aircraft_cases{end+1} = make_aircraft_from_matrix('my-B', B, umin, umax, '');
% If cpp_tag is a formal case name, test.m loads results/cpp_outputs/output_cpp_<cpp_tag>_*.csv for C++ comparison.
% If cpp_tag is '', only MATLAB outputs are simulated and saved/plotted.
aircraft_cases = select_aircraft_cases(aircraft_cases, aircraft_case_selection);

%% Tie-break settings shared with alloc_cpp/src/ControlAllocation/ControlAllocation.h
tie_opts.tie_rel_tol = 1e-5;
tie_opts.tie_abs_tol = 1e-6;
tie_opts.zero_tie_abs_tol = 3e-5;

fprintf('MATLAB PCA methods: pca_dir=PCA/DP_LPCA, pca_dpscaled=PCA/DPscaled_LPCA, pca_prio=PCA/DP_LPCA_prio + restoring_cpp\n');
fprintf('C++ methods: cpp_dir=DP_LPCA, cpp_dpscaled=DPscaled_LPCA, cpp_prio=DP_LPCA_prio + restoring() aligned to restoring_cpp\n');
fprintf('Restoring: %s\n', string(use_restoring));
fprintf('Allocator methods enabled: %s\n', strjoin(allocation_methods_to_run, ', '));
fprintf('Aircraft cases enabled: %d\n', numel(aircraft_cases));
for case_idx = 1:numel(aircraft_cases)
    print_aircraft_case_model(aircraft_cases{case_idx});
end

%% Simulate flight process
tic;
results = cell(size(aircraft_cases));
for case_idx = 1:numel(aircraft_cases)
    results{case_idx} = simulate_flight_process(aircraft_cases{case_idx}, v_test, tie_opts, allocation_methods_to_run, use_restoring);
end
elapsed_time = toc;
fprintf('MATLAB reference execution time: %.2f 秒\n', elapsed_time);
fprintf('Test sample window: %.6g to %.6g s (%d samples)\n', t(1), t(end), len_command_px4);

%% C++ outputs
if any(cellfun(@(result) isfield(result, 'cpp_tag') && ~isempty(result.cpp_tag), results))
    ensure_cpp_outputs();
end

for case_idx = 1:numel(results)
    results{case_idx} = load_cpp_outputs_for_case(results{case_idx}, test_sample_indices, use_restoring);
end

%% Reports
for case_idx = 1:numel(results)
    report_case(results{case_idx}, command_px4, len_command_px4, t);
    if has_report(reports_to_run, 'allocator_summary')
        report_allocator_methods(results{case_idx}, allocation_methods_to_run, command_px4, t, ...
                                 allocator_sample_print_count > 0, allocator_sample_print_count);
    end

    if has_report(reports_to_run, 'method_diff')
        report_allocator_method_differences(results{case_idx}, allocation_methods_to_run);
    end
end

if has_report(reports_to_run, 'case_diff')
    report_case_comparisons(results, case_comparison_pairs, allocation_methods_to_run, t);
end

if has_report(reports_to_run, 'reachable_set')
    plot_reachable_sets_for_results(results);
end

%% Save once, then rerun only plot_test_results.m when changing plots
result4 = [];
result6 = [];
if numel(results) >= 1
    result4 = results{1};
end
if numel(results) >= 2
    result6 = results{2};
end
save('test_results.mat', 'results', 'result4', 'result6', 'command_px4', 'len_command_px4', 't', 'tie_opts', 'test_time_window_s', 'test_sample_indices', ...
     'allocation_method_catalog', 'allocation_method_selection', 'allocation_methods_to_run', 'aircraft_case_selection', 'case_comparison_pairs', ...
     'reports_to_run', 'allocator_sample_print_count', 'use_restoring');
plot_test_results('test_results.mat');

function aircraft = make_aircraft_4()
    % Four-effector ducted-fan model.  This form keeps the geometry and
    % inertia parameters visible; it is numerically equal to alloc_cpp/test/main.cpp.
    %
    % 行满秩时 pinv(B) = The Moore-Penrose Pseudo-inverse
    % If B is full row rank, The Moore-Penrose Pseudo-inverse
    % B^+ = B^T (B B^T)^{-1}, since
    % B = K*P,
    % K = I\diag([l1 l1 l2])*k_omega2force,
    % P = [-1     0     1     0;
    %       0    -1     0     1;
    %       1     1     1     1];
    % K = diag([ k_omega2force*l1/I_x  k_omega2force*l1/I_y  k_omega2force*l2/I_z ])
    %
    % B^{dagger} = P^T K K^{-1} (P P^T)^{-1} K^{-1}
    %            = P^T (P P^T)^{-1} K^{-1}
    %            = P^{dagger} K^{-1}
    %
    % P^{dagger} = [-0.5000  -0.0000   0.2500;
    %                0       -0.5000   0.2500;
    %                0.5000  -0.0000   0.2500;
    %                0        0.5000   0.2500]
    %
    % B^{dagger} = P^{dagger} K^{-1}
    %            = P^{dagger} * diag([ I_x/(k_omega2force*l1)  I_y/(k_omega2force*l1)  I_z/(k_omega2force*l2) ])

    l1 = 0.167;
    l2 = 0.069;
    k_omega2force = 1;
    I_x = 0.01149;
    I_y = 0.01153;
    I_z = 0.00487;
    I = diag([I_x, I_y, I_z]);

    %=============================4==================================
    B = I \ [-l1     0       l1     0;
              0      -l1     0      l1;
              l2     l2      l2     l2] * k_omega2force;
    
    % B = I\diag([l1 l1 l2])[-1     0     1     0;
    %                         0    -1     0     1;
    %                         1     1     1     1]*k_omega2force;
    %
    % B^+ = B^T (B B^T)^{-1}
    %    1X
    % 4     2Y
    %    3

    ulim = 0.3491;
    aircraft.name = '15003_ductedfan4';
    aircraft.cpp_tag = aircraft.name;
    aircraft.B = B;
    aircraft.umin = ones(4, 1) * -ulim;
    aircraft.umax = ones(4, 1) * ulim;
end

function aircraft = make_aircraft_6()
    % SHC09 six-effector model, same physical construction as alloc_cpp/test/main.cpp.
    l1 = 0.292166;
    l2 = 0.073699;
    k_omega2force = 1.93;
    I_x = 0.0438;
    I_y = 0.0436;
    I_z = 0.005006;
    d = 60*pi/180;
    I = diag([I_x, I_y, I_z]);

    %=============================6==================================
    % Original geometry note:
    % B = I\[-l1 -l1*cos(d) l1*cos(d) l1 l1*cos(d) -l1*cos(d);
    %        0 -l1*sin(d) -l1*sin(d) 0 l1*sin(d) l1*sin(d);
    %        l2 l2 l2 l2 l2 l2]*k_v;
    %      1X
    %  6       2
    % -------------->Y
    %  5       3
    %      4
    %
    % The active sign convention below matches alloc_cpp/test/main.cpp.
    B = I \ [-l1, -l1*cos(d),  l1*cos(d),  l1,  l1*cos(d), -l1*cos(d);
              0,   l1*sin(d),  l1*sin(d),  0,  -l1*sin(d), -l1*sin(d);
              l2,  l2,         l2,         l2,  l2,         l2] * k_omega2force;
    % K = I\diag([l1 l1 l2])*k_omega2force;
    % P = [-1     -cos(d)     cos(d)    1    cos(d)  -cos(d);
    %       0,    sin(d),     sin(d),   0,  -sin(d), -sin(d);
    %       1     1           1         1    1        1];%B = K*P
    % B_inv_pid=[-1 0 1;-1 1 1;1 1 1;1 0 1;1 -1 1;-1 -1 1]
    % B_pid=pinv(B_inv_pid)=[-0.1667   -0.1667    0.1667    0.1667    0.1667   -0.1667;
    %                         0         0.2500    0.2500    0        -0.2500   -0.2500;
    %                         0.1667    0.1667    0.1667    0.1667    0.1667    0.1667];
    ulim = 40*pi/180;
    aircraft.name = '15006_SHC09';
    aircraft.cpp_tag = aircraft.name;
    aircraft.B = B;
    aircraft.umin = ones(6, 1) * -ulim;
    aircraft.umax = ones(6, 1) * ulim;
end

function aircraft = make_aircraft_from_matrix(name, B, umin, umax, cpp_tag)
    if nargin < 5
        cpp_tag = '';
    end
    aircraft.name = name;
    aircraft.cpp_tag = cpp_tag;
    aircraft.B = B;
    aircraft.umin = umin;
    aircraft.umax = umax;
end

function catalog = get_allocation_method_catalog()
    % Method name catalog. Keep this list in one place so the user-facing
    % selection above can be 'all', numeric indices, or names.
    catalog = { ...
        'inv', ...                    % pinv(B)*y + clamp + restoring_cpp
        'pca_dir', ...                % /PCA/DP_LPCA + restoring_cpp
        'pca_dpscaled', ...           % /PCA/DPscaled_LPCA + restoring_cpp
        'pca_prio', ...               % /PCA/DP_LPCA_prio + restoring_cpp
        'cpp_dir', ...                % alloc_cpp/test/main.cpp DP_LPCA CSV output
        'cpp_dpscaled', ...           % alloc_cpp/test/main.cpp DPscaled_LPCA CSV output
        'cpp_prio', ...               % alloc_cpp/test/main.cpp DP_LPCA_prio CSV output
        'lib_lpwrap_db', ...          % control_allocation_lib LPwrap.m LPmethod=0
        'lib_lpwrap_dbinf', ...       % control_allocation_lib LPwrap.m LPmethod=1
        'lib_lpwrap_dir', ...          % control_allocation_lib LPwrap.m LPmethod=2
        'lib_lpwrap_dpscaled', ...    % control_allocation_lib LPwrap.m LPmethod=3
        'lib_lpwrap_mo', ...          % control_allocation_lib LPwrap.m LPmethod=4
        'lib_lpwrap_sb', ...          % control_allocation_lib LPwrap.m LPmethod=5
        'lib_lpwrap_par_db', ...      % control_allocation_lib LPwrap_par.m LPmethod=0
        'lib_lpwrap_par_dbinf', ...   % control_allocation_lib LPwrap_par.m LPmethod=1
        'lib_lpwrap_par_dir', ...      % control_allocation_lib LPwrap_par.m LPmethod=2
        'lib_lpwrap_par_dpscaled', ...% control_allocation_lib LPwrap_par.m LPmethod=3
        'lib_lpwrap_par_mo', ...      % control_allocation_lib LPwrap_par.m LPmethod=4
        'lib_lpwrap_par_sb', ...      % control_allocation_lib LPwrap_par.m LPmethod=5
        'lib_lpwrap_incre', ...       % control_allocation_lib LPwrap.m LPmethod=2, incremental-form placeholder
        'lib_cgiwrap', ...            % control_allocation_lib CGIwrap.m
        'lib_dawrap', ...             % control_allocation_lib DAwrap.m
        'lib_vjawrap', ...            % control_allocation_lib VJAwrap.m
        'wls', ...                    % QCAT weighted least-squares allocator
        'wls_gen', ...                % Generated WLS allocator
        'dir_linprog', ...            % Direct allocation using linprog
        'dir_linprog_re', ...         % Reformulated direct allocation using linprog
        'dir_linprog_re_bound', ...   % Reformulated bounded direct allocation
        'use_lp_lib', ...             % Standard-form LP-library experiment
        'allocator_dir_lpwrap_4'};    % 4-effector generated wrapper only
end

function methods = resolve_allocation_method_selection(selection, catalog)
    if nargin < 2 || isempty(catalog)
        catalog = {'inv', 'pca_dir', 'pca_dpscaled', 'pca_prio', 'wls'};
    end

    if isempty(selection)
        methods = validate_allocation_methods(catalog, catalog);
        return;
    end

    if ischar(selection) || isstring(selection)
        selection = cellstr(selection);

        if numel(selection) == 1 && strcmpi(selection{1}, 'all')
            methods = validate_allocation_methods(catalog, catalog);
            return;
        end

        methods = validate_allocation_methods(selection, catalog);
        return;
    end

    if isnumeric(selection)
        methods = validate_allocation_methods(catalog(selection), catalog);
        return;
    end

    if iscell(selection)
        methods = validate_allocation_methods(selection, catalog);
        return;
    end

    error('test:BadAllocatorSelection', 'allocation_method_selection must be ''all'', numeric indices, or method names.');
end

function methods = validate_allocation_methods(methods, catalog)
    methods = cellfun(@(method) lower(char(method)), methods, 'UniformOutput', false);
    unknown = setdiff(methods, catalog);
    if ~isempty(unknown)
        error('test:UnknownAllocatorMethod', ...
              'Unknown allocator method(s): %s. Use one of: %s', ...
              strjoin(unknown, ', '), strjoin(catalog, ', '));
    end
    methods = unique(methods, 'stable');
end

function selected_cases = select_aircraft_cases(aircraft_cases, selection)
    if isempty(selection)
        selected_cases = aircraft_cases;
        return;
    end

    if ischar(selection) || isstring(selection)
        selection = cellstr(selection);

        if numel(selection) == 1 && strcmpi(selection{1}, 'all')
            selected_cases = aircraft_cases;
            return;
        end
    end

    if isnumeric(selection)
        selected_cases = aircraft_cases(selection);
        return;
    end

    if ~iscell(selection)
        error('test:BadCaseSelection', 'aircraft_case_selection must be ''all'', numeric indices, or case names.');
    end

    selected_cases = cell(1, numel(selection));

    for selection_idx = 1:numel(selection)
        wanted_name = selection{selection_idx};
        found_idx = -1;

        for case_idx = 1:numel(aircraft_cases)
            if strcmp(aircraft_cases{case_idx}.name, wanted_name)
                found_idx = case_idx;
                break;
            end
        end

        if found_idx < 1
            available_names = cellfun(@(aircraft) aircraft.name, aircraft_cases, 'UniformOutput', false);
            error('test:UnknownAircraftCase', 'Unknown aircraft case "%s". Available cases: %s', ...
                  wanted_name, strjoin(available_names, ', '));
        end

        selected_cases{selection_idx} = aircraft_cases{found_idx};
    end
end

function enabled = has_report(reports_to_run, report_name)
    if ischar(reports_to_run) || isstring(reports_to_run)
        reports_to_run = cellstr(reports_to_run);
    end

    enabled = any(strcmpi(reports_to_run, 'all')) || any(strcmpi(reports_to_run, report_name));
end

function print_aircraft_case_model(aircraft)
    fprintf('\n=== %s model ===\n', aircraft.name);
    fprintf('B =\n');
    disp(aircraft.B);
    fprintf('umin = [');
    fprintf(' %.6f', aircraft.umin);
    fprintf(' ]^T\n');
    fprintf('umax = [');
    fprintf(' %.6f', aircraft.umax);
    fprintf(' ]^T\n');
end

function aircraft = make_shw09_vtol_mc_full8()
    model = make_shw09_vtol_mc_model();

    aircraft.name = '15008_SHW09_vtol_full8';
    aircraft.cpp_tag = '';
    aircraft.B = model.B;
    aircraft.umin = model.umin;
    aircraft.umax = model.umax;
    aircraft.plot_num_outputs = 8;
    aircraft.plot_active_idx = 1:8;
end

function aircraft = make_shw09_vtol_mc_reduced6()
    model = make_shw09_vtol_mc_model();
    [B_reduced, umin_reduced, umax_reduced, active_idx] = active_reachable_view_model(model.B, model.umin, model.umax);

    aircraft.name = '15008_SHW09_vtol_reduced6';
    aircraft.cpp_tag = '';
    aircraft.B = B_reduced;
    aircraft.umin = umin_reduced;
    aircraft.umax = umax_reduced;
    aircraft.plot_num_outputs = 8;
    aircraft.plot_active_idx = active_idx;
end

function model = make_shw09_vtol_mc_model()
    I_x = 0.050636;
    I_y = 0.042954;
    I_z = 0.012668;
    L_1 = 0.42;
    L_2 = 0.073699;
    k = 4.2;
    d = 60*pi/180;
    ulim = 0.6981;

    B = zeros(3, 8);
    B(1,1:6) = [-L_1, -cos(d)*L_1, cos(d)*L_1, L_1, cos(d)*L_1, -cos(d)*L_1] * k / I_x;
    B(2,1:6) = [0, sin(d)*L_1, sin(d)*L_1, 0, -sin(d)*L_1, -sin(d)*L_1] * k / I_y;
    B(3,1:6) = ones(1,6) * L_2 * k / I_z;

    umin = ones(8,1) * -ulim;
    umax = ones(8,1) *  ulim;
    umin(7:8) = 0;
    umax(7:8) = 0;

    model.B = B;
    model.umin = umin;
    model.umax = umax;
end

function result = simulate_flight_process(aircraft, v, tie_opts, allocation_methods_to_run, use_restoring)
    if nargin < 4 || isempty(allocation_methods_to_run)
        allocation_methods_to_run = {'pca_dir', 'pca_dpscaled', 'pca_prio'};
    end
    if nargin < 5
        use_restoring = true;
    end

    B = aircraft.B;
    umin = aircraft.umin;
    umax = aircraft.umax;
    [k, m] = size(B);
    [~, N] = size(v);

    result = aircraft;
    result.allocation_methods = allocation_methods_to_run;
    result.use_restoring = use_restoring;

    for method_idx = 1:numel(allocation_methods_to_run)
        method = allocation_methods_to_run{method_idx};
        result.alloc.(method).u_raw = nan(m, N);
        result.alloc.(method).u = nan(m, N);
        result.alloc.(method).y_achieved = nan(k, N);
        result.alloc.(method).ok = false(1, N);
        result.alloc.(method).err_msg = strings(1, N);
    end

    %% setup function of allocation lib
    % LPwrap/CGIwrap/DAwrap/VJAwrap use this global to size fallback zero outputs.
    global NumU
    NumU=m;

    %% simulate flight process
    for idx = 1:N  % or x:N for debug
        command = v(:, idx);

        for method_idx = 1:numel(allocation_methods_to_run)
            method = allocation_methods_to_run{method_idx};

            if is_cpp_allocator_method(method)
                % C++ allocator outputs are loaded from results/cpp_outputs
                % after MATLAB simulation, so they are registered into
                % result.alloc by load_cpp_outputs_for_case().
                continue;
            end

            [u_raw, u_alloc, ok, err_msg] = run_allocator_method(method, command, B, umin, umax, tie_opts, use_restoring);

            result.([method '_raw'])(:, idx) = u_raw;
            result.(method)(:, idx) = u_alloc;
            result.alloc.(method).u_raw(:, idx) = u_raw;
            result.alloc.(method).u(:, idx) = u_alloc;
            result.alloc.(method).ok(idx) = ok;
            result.alloc.(method).err_msg(idx) = err_msg;

            if ok
                result.alloc.(method).y_achieved(:, idx) = B * u_alloc;
            end

            switch method
                case 'inv'
                    result.inv_raw(:, idx) = u_raw;
                    result.inv(:, idx) = u_alloc;

                case 'wls'
                    result.wls_raw(:, idx) = u_raw;
                    result.wls(:, idx) = u_alloc;
            end
        end

    end
end

function tf = is_cpp_allocator_method(method)
    tf = strncmpi(char(method), 'cpp_', 4);
end

function [u_raw, u, ok, err_msg] = run_allocator_method(method, y, B, umin, umax, tie_opts, use_restoring)
    if nargin < 7
        use_restoring = true;
    end
    [k, m] = size(B);
    u_raw = nan(m, 1);
    u = nan(m, 1);
    ok = true;
    err_msg = "";

    try
        switch lower(method)
            case 'inv'
                u_raw = clamp_u(pinv(B) * y, umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'pca_dir'
                [u_tmp, ~, ~] = DP_LPCA(y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'pca_dpscaled'
                [u_tmp, ~, ~, ~] = DPscaled_LPCA(y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'pca_prio'
                m_higher = zeros(k, 1);
                [u_tmp, ~, ~] = DP_LPCA_prio(m_higher, y, B, umin, umax, 100, tie_opts);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_db'
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 0));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_dbinf'
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 1));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_dir'
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 2));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_dpscaled'
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 3));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_mo'
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 4));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_sb'
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin, umax, 5));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_par_db'
                u_tmp = LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 0), y, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_par_dbinf'
                u_tmp = LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 1), y, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_par_dir'
                u_tmp = LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 2), y, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_par_dpscaled'
                u_tmp = LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 3), y, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_par_mo'
                u_tmp = LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 4), y, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_par_sb'
                u_tmp = LPwrap_par(make_lpwrap_in_mat(B, y, umin, umax, 5), y, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_lpwrap_incre'
                u0 = zeros(m, 1);
                u_tmp = LPwrap(make_lpwrap_in_mat(B, y, umin - u0, umax - u0, 2));
                u_raw = clamp_u(u_tmp(:), umin - u0, umax - u0) + u0;
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_cgiwrap'
                u_tmp = CGIwrap(make_book_wrapper_in_mat(B, y, umin, umax));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_dawrap'
                u_tmp = DAwrap(make_book_wrapper_in_mat(B, y, umin, umax));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'lib_vjawrap'
                u_tmp = VJAwrap(make_book_wrapper_in_mat(B, y, umin, umax));
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'wls'
                Wv = eye(k);
                Wu = eye(m);
                ud = zeros(m, 1);
                gam = 1e6;
                u0 = (umin + umax) / 2;
                W0 = zeros(m, 1);
                [u_tmp, ~, ~] = wls_alloc(B, y, umin, umax, Wv, Wu, ud, gam, u0, W0, 100);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'wls_gen'
                Wv = eye(k);
                Wu = eye(m);
                ud = zeros(m, 1);
                gam = 1e6;
                u0 = zeros(m, 1);
                W0 = zeros(m, 1);
                u_tmp = wls_alloc_gen(B, y, umin, umax, Wv, Wu, ud, gam, u0, W0, 100, m);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'dir_linprog'
                [u_tmp, ~] = dir_alloc_linprog(B, y, umin, umax, 1e4);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'dir_linprog_re'
                [u_tmp, ~] = dir_alloc_linprog_re(B, y, umin, umax);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'dir_linprog_re_bound'
                [u_tmp, ~] = dir_alloc_linprog_re_bound(B, y, umin, umax, 1e4);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'use_lp_lib'
                [u_tmp, ~] = use_LP_lib(B, y, umin, umax);
                u_raw = clamp_u(u_tmp(:), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            case 'allocator_dir_lpwrap_4'
                if m ~= 4
                    error('test:AllocatorDimensionMismatch', 'allocator_dir_LPwrap_4 only supports m=4, current m=%d.', m);
                end
                [u_tmp, ~, ~] = allocator_dir_LPwrap_4(single(B), single(y), single(umin), single(umax));
                u_raw = clamp_u(double(u_tmp(:)), umin, umax);
                u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring);

            otherwise
                error('test:UnknownAllocatorMethod', 'Unknown allocator method: %s', method);
        end
    catch ME
        ok = false;
        err_msg = string(ME.message);
    end
end

function IN_MAT = make_lpwrap_in_mat(B, y, umin, umax, lp_method)
    [~, m] = size(B);
    active_effectors = ones(1, m);
    IN_MAT = [B                y(:);
              umin(:)'         0;
              umax(:)'         0;
              active_effectors lp_method];
end

function IN_MAT = make_book_wrapper_in_mat(B, y, umin, umax)
    % CGIwrap/DAwrap/VJAwrap use the same [B y; umin 0; umax 0; INDX 0]
    % matrix shape as LPwrap, but the last scalar is not an LPmethod.
    IN_MAT = make_lpwrap_in_mat(B, y, umin, umax, 0);
end

function u = clamp_u(u, umin, umax)
    u = min(max(u(:), umin), umax);
end

function u = apply_optional_restoring(B, u_raw, umin, umax, use_restoring)
    if should_apply_restoring(use_restoring)
        u = restoring_cpp(B, u_raw, umin, umax);
    else
        u = u_raw;
    end
end

function tf = should_apply_restoring(use_restoring)
    if ischar(use_restoring) || isstring(use_restoring)
        mode = lower(char(use_restoring));
        tf = any(strcmp(mode, {'true', 'on', 'yes', 'restoring', 'both'}));
    else
        tf = logical(use_restoring);
    end
end

function ensure_cpp_outputs()
    cpp_output_dir = fullfile('results', 'cpp_outputs');
    cpp_output_files = [
        cpp_output_file_names(cpp_output_dir, '15003_ductedfan4'), ...
        cpp_output_file_names(cpp_output_dir, '15006_SHC09')];
    if ~isfile('./alloc_cpp/build/main')
        error('test:MissingCppBinary', 'Build alloc_cpp/build/main before running test.m.');
    end
    if all(isfile(cpp_output_files))
        cpp_binary_info = dir('./alloc_cpp/build/main');
        output_info = dir(cpp_output_files(1));
        if ~isempty(cpp_binary_info) && ~isempty(output_info) && output_info.datenum >= cpp_binary_info.datenum
            return;
        end
    end
    status = system('./alloc_cpp/build/main');
    if status ~= 0
        error('test:CppRunFailed', 'Failed to run alloc_cpp/build/main.');
    end
end

function files = cpp_output_file_names(cpp_output_dir, cpp_tag)
    files = [fullfile(cpp_output_dir, sprintf('output_cpp_%s_dir.csv', cpp_tag)), ...
             fullfile(cpp_output_dir, sprintf('output_cpp_%s_dir_raw.csv', cpp_tag)), ...
             fullfile(cpp_output_dir, sprintf('output_cpp_%s_dpscaled.csv', cpp_tag)), ...
             fullfile(cpp_output_dir, sprintf('output_cpp_%s_dpscaled_raw.csv', cpp_tag)), ...
             fullfile(cpp_output_dir, sprintf('output_cpp_%s_prio.csv', cpp_tag)), ...
             fullfile(cpp_output_dir, sprintf('output_cpp_%s_prio_raw.csv', cpp_tag))];
end

function result = load_cpp_outputs_for_case(result, sample_indices, use_restoring)
    if nargin < 3
        use_restoring = true;
    end
    if ~isfield(result, 'cpp_tag') || isempty(result.cpp_tag)
        return;
    end

    tag = result.cpp_tag;
    cpp_output_dir = fullfile('results', 'cpp_outputs');
    dir_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_dir.csv', tag));
    dir_raw_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_dir_raw.csv', tag));
    dpscaled_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_dpscaled.csv', tag));
    dpscaled_raw_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_dpscaled_raw.csv', tag));
    prio_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_prio.csv', tag));
    prio_raw_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_prio_raw.csv', tag));

    if isfile(dir_file)
        result.cpp_dir = slice_cpp_output(readmatrix(dir_file)', sample_indices);
    end
    if isfile(dir_raw_file)
        result.cpp_dir_raw = slice_cpp_output(readmatrix(dir_raw_file)', sample_indices);
    end
    if isfile(dpscaled_file)
        result.cpp_dpscaled = slice_cpp_output(readmatrix(dpscaled_file)', sample_indices);
    end
    if isfile(dpscaled_raw_file)
        result.cpp_dpscaled_raw = slice_cpp_output(readmatrix(dpscaled_raw_file)', sample_indices);
    end
    if isfile(prio_file)
        result.cpp_prio = slice_cpp_output(readmatrix(prio_file)', sample_indices);
    end
    if isfile(prio_raw_file)
        result.cpp_prio_raw = slice_cpp_output(readmatrix(prio_raw_file)', sample_indices);
    end

    if ~should_apply_restoring(use_restoring)
        if isfield(result, 'cpp_dir_raw')
            result.cpp_dir = result.cpp_dir_raw;
        end
        if isfield(result, 'cpp_dpscaled_raw')
            result.cpp_dpscaled = result.cpp_dpscaled_raw;
        end
        if isfield(result, 'cpp_prio_raw')
            result.cpp_prio = result.cpp_prio_raw;
        end
    end

    if isfield(result, 'cpp_dir_raw') && isfield(result, 'cpp_dir')
        result = attach_cpp_allocator_output(result, 'cpp_dir', result.cpp_dir_raw, result.cpp_dir);
    end
    if isfield(result, 'cpp_dpscaled_raw') && isfield(result, 'cpp_dpscaled')
        result = attach_cpp_allocator_output(result, 'cpp_dpscaled', result.cpp_dpscaled_raw, result.cpp_dpscaled);
    end
    if isfield(result, 'cpp_prio_raw') && isfield(result, 'cpp_prio')
        result = attach_cpp_allocator_output(result, 'cpp_prio', result.cpp_prio_raw, result.cpp_prio);
    end
end

function x = slice_cpp_output(x_full, sample_indices)
    if max(sample_indices) <= size(x_full, 2)
        x = x_full(:, sample_indices);
    else
        % If the C++ binary was explicitly run on an already-sliced input,
        % keep that output rather than failing here.
        x = x_full;
    end
end

function result = attach_cpp_allocator_output(result, method, u_raw, u)
    result.([method '_raw']) = u_raw;
    result.(method) = u;

    result.alloc.(method).u_raw = u_raw;
    result.alloc.(method).u = u;
    result.alloc.(method).y_achieved = result.B * u;
    result.alloc.(method).ok = all(isfinite(u), 1);
    result.alloc.(method).err_msg = strings(1, size(u, 2));
    result.alloc.(method).err_msg(~result.alloc.(method).ok) = "non-finite C++ output";
end

function report_case(result, command_px4, len_command_px4, t)
    method_names = {};
    if isfield(result, 'alloc')
        method_names = fieldnames(result.alloc);
    end
    final_label = 'restored';
    if isfield(result, 'use_restoring') && ~should_apply_restoring(result.use_restoring)
        final_label = 'no restoring';
    end
    show_raw_report = ~isfield(result, 'use_restoring') || should_apply_restoring(result.use_restoring);

    if isempty(method_names) || any(strcmp(method_names, 'pca_dir'))
        if show_raw_report && isfield(result, 'cpp_dir_raw')
            report_cpp_matlab_case([result.name ' pca_dir raw'], result.B, result.cpp_dir_raw, result.pca_dir_raw, command_px4, len_command_px4, t);
        end
        if isfield(result, 'cpp_dir')
            report_cpp_matlab_case([result.name ' pca_dir ' final_label], result.B, result.cpp_dir, result.pca_dir, command_px4, len_command_px4, t);
        else
            report_matlab_case([result.name ' pca_dir ' final_label], result.B, result.pca_dir, command_px4, len_command_px4);
        end
    end

    if isempty(method_names) || any(strcmp(method_names, 'pca_dpscaled'))
        if show_raw_report && isfield(result, 'cpp_dpscaled_raw')
            report_cpp_matlab_case([result.name ' pca_dpscaled raw'], result.B, result.cpp_dpscaled_raw, result.pca_dpscaled_raw, command_px4, len_command_px4, t);
        end
        if isfield(result, 'cpp_dpscaled')
            report_cpp_matlab_case([result.name ' pca_dpscaled ' final_label], result.B, result.cpp_dpscaled, result.pca_dpscaled, command_px4, len_command_px4, t);
        else
            report_matlab_case([result.name ' pca_dpscaled ' final_label], result.B, result.pca_dpscaled, command_px4, len_command_px4);
        end
    end

    if isempty(method_names) || any(strcmp(method_names, 'pca_prio'))
        if show_raw_report && isfield(result, 'cpp_prio_raw')
            report_cpp_matlab_case([result.name ' pca_prio raw'], result.B, result.cpp_prio_raw, result.pca_prio_raw, command_px4, len_command_px4, t);
        end
        if isfield(result, 'cpp_prio')
            report_cpp_matlab_case([result.name ' pca_prio ' final_label], result.B, result.cpp_prio, result.pca_prio, command_px4, len_command_px4, t);
        else
            report_matlab_case([result.name ' pca_prio ' final_label], result.B, result.pca_prio, command_px4, len_command_px4);
        end
    end

end

function report_cpp_matlab_case(name, B, x_cpp_full, x_matlab_full, command_full, len_command_px4, t)
    n = min([len_command_px4, size(x_cpp_full, 2), size(x_matlab_full, 2), size(command_full, 2)]);
    x_cpp = x_cpp_full(:, 1:n);
    x_matlab = x_matlab_full(:, 1:n);
    command = command_full(:, 1:n);

    U_cpp = B * x_cpp;
    U_matlab = B * x_matlab;

    du = max(abs(x_cpp - x_matlab), [], 'all');
    dU = max(abs(U_cpp - U_matlab), [], 'all');
    e_cpp = max(abs(U_cpp - command), [], 'all');
    e_matlab = max(abs(U_matlab - command), [], 'all');

    fprintf('\n[%s]\n', name);
    fprintf('max |u_cpp - u_matlab| = %.6g\n', du);
    fprintf('max |B*u_cpp - B*u_matlab| = %.6g\n', dU);
    fprintf('max |B*u_cpp - command| = %.6g\n', e_cpp);
    fprintf('max |B*u_matlab - command| = %.6g\n', e_matlab);

    [~, worst_linear_idx] = max(abs(x_cpp - x_matlab), [], 'all', 'linear');
    [worst_u_idx, worst_t_idx] = ind2sub(size(x_cpp), worst_linear_idx);
    fprintf('worst u index = %d, sample = %d, time = %.6g s\n', worst_u_idx, worst_t_idx, t(worst_t_idx));
    fprintf('command(:,sample) = ['); fprintf(' %.12g', command(:, worst_t_idx)); fprintf(' ]^T\n');
    fprintf('u_cpp(:,sample)   = ['); fprintf(' %.12g', x_cpp(:, worst_t_idx)); fprintf(' ]^T\n');
    fprintf('u_matlab(:,sample)= ['); fprintf(' %.12g', x_matlab(:, worst_t_idx)); fprintf(' ]^T\n');
end

function report_matlab_case(name, B, x_matlab_full, command_full, len_command_px4)
    n = min([len_command_px4, size(x_matlab_full, 2), size(command_full, 2)]);
    x_matlab = x_matlab_full(:, 1:n);
    command = command_full(:, 1:n);
    U_matlab = B * x_matlab;
    e_matlab = max(abs(U_matlab - command), [], 'all');

    fprintf('\n[%s]\n', name);
    fprintf('max |B*u_matlab - command| = %.6g\n', e_matlab);
end

function report_allocator_methods(result, methods_to_run, command_full, t, print_samples, sample_count)
    if ~isfield(result, 'alloc')
        return;
    end

    n = min([size(command_full, 2), numel(t)]);
    fprintf('\n=== MATLAB allocator smoke: %s ===\n', result.name);

    active_limits = abs(result.umax - result.umin) > 1e-9;
    if ~any(active_limits)
        active_limits = true(size(result.umin));
    end

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};

        if ~isfield(result.alloc, method)
            fprintf('[%-8s] missing\n', method);
            continue;
        end

        data = result.alloc.(method);
        ok = data.ok(1:n);
        fail_count = n - nnz(ok);

        if any(ok)
            y = data.y_achieved(:, 1:n);
            command = command_full(:, 1:n);
            err = y(:, ok) - command(:, ok);
            err_norm = vecnorm(err, 2, 1);
            u = data.u(:, 1:n);
            u_ok = u(:, ok);
            sat_all = mean(any(abs(u_ok - result.umin) < 1e-6 | abs(u_ok - result.umax) < 1e-6, 1)) * 100;
            sat_active = mean(any(abs(u_ok(active_limits,:) - result.umin(active_limits)) < 1e-6 | ...
                                  abs(u_ok(active_limits,:) - result.umax(active_limits)) < 1e-6, 1)) * 100;

            fprintf('[%-8s] ok=%d/%d, fail=%d, max|err|=%.6g, rms|err|=%.6g, max|u|=%.6g, sat_active=%.2f%%, sat_all=%.2f%%\n', ...
                    method, nnz(ok), n, fail_count, max(abs(err), [], 'all'), sqrt(mean(err_norm.^2)), max(abs(u_ok), [], 'all'), sat_active, sat_all);
        else
            fprintf('[%-8s] ok=0/%d, fail=%d\n', method, n, fail_count);
        end

        if fail_count > 0
            first_fail = find(~ok, 1);
            fprintf('          first failure sample=%d, t=%.6g s, error=%s\n', first_fail, t(first_fail), data.err_msg(first_fail));
        end

        if print_samples
            print_count = min(sample_count, n);
            for sample_idx = 1:print_count
                fprintf('          sample=%d, t=%.6g s, ok=%d\n', sample_idx, t(sample_idx), data.ok(sample_idx));
                fprintf('            y  = ['); fprintf(' %.9g', command_full(:, sample_idx)); fprintf(' ]^T\n');
                fprintf('            u  = ['); fprintf(' %.9g', data.u(:, sample_idx)); fprintf(' ]^T\n');
                fprintf('            Bu = ['); fprintf(' %.9g', data.y_achieved(:, sample_idx)); fprintf(' ]^T\n');
            end
        end
    end
end

function report_allocator_method_differences(result, methods_to_run)
    if ~isfield(result, 'alloc')
        return;
    end

    base_method = '';
    if isfield(result.alloc, 'inv') && any(result.alloc.inv.ok)
        base_method = 'inv';
    else
        for method_idx = 1:numel(methods_to_run)
            method = methods_to_run{method_idx};
            if isfield(result.alloc, method) && any(result.alloc.(method).ok)
                base_method = method;
                break;
            end
        end
    end

    if isempty(base_method)
        fprintf('\n--- method differences in %s ---\n', result.name);
        fprintf('No successful baseline allocator result.\n');
        return;
    end

    fprintf('\n--- method differences in %s, baseline=%s ---\n', result.name, base_method);
    base = result.alloc.(base_method);

    for method_idx = 1:numel(methods_to_run)
        method = methods_to_run{method_idx};
        if strcmp(method, base_method) || ~isfield(result.alloc, method)
            continue;
        end

        data = result.alloc.(method);
        n = min([numel(base.ok), numel(data.ok)]);
        common_ok = base.ok(1:n) & data.ok(1:n);

        if ~any(common_ok)
            fprintf('[%-8s vs %-8s] no common successful samples\n', method, base_method);
            continue;
        end

        du = max(abs(data.u(:, common_ok) - base.u(:, common_ok)), [], 'all');
        dy = max(abs(data.y_achieved(:, common_ok) - base.y_achieved(:, common_ok)), [], 'all');
        fprintf('[%-8s vs %-8s] common=%d, max|du|=%.6g, max|d(Bu)|=%.6g\n', ...
                method, base_method, nnz(common_ok), du, dy);
    end
end

function report_case_comparisons(results, case_comparison_pairs, methods_to_run, t)
    if isempty(case_comparison_pairs)
        return;
    end

    for pair_idx = 1:size(case_comparison_pairs, 1)
        name_a = case_comparison_pairs{pair_idx, 1};
        name_b = case_comparison_pairs{pair_idx, 2};
        idx_a = find_result_index_by_name(results, name_a);
        idx_b = find_result_index_by_name(results, name_b);

        fprintf('\n=== case comparison: %s vs %s ===\n', name_a, name_b);

        if idx_a < 1 || idx_b < 1
            fprintf('Skipped: one or both cases are not enabled.\n');
            continue;
        end

        result_a = results{idx_a};
        result_b = results{idx_b};
        n_act = max(get_plot_num_outputs(result_a), get_plot_num_outputs(result_b));

        for method_idx = 1:numel(methods_to_run)
            method = methods_to_run{method_idx};

            if ~isfield(result_a.alloc, method) || ~isfield(result_b.alloc, method)
                fprintf('[%-8s] missing in one case\n', method);
                continue;
            end

            data_a = result_a.alloc.(method);
            data_b = result_b.alloc.(method);
            n = min([numel(data_a.ok), numel(data_b.ok), numel(t)]);
            common_ok = data_a.ok(1:n) & data_b.ok(1:n);

            if ~any(common_ok)
                fprintf('[%-8s] no common successful samples\n', method);
                continue;
            end

            u_a = embed_u_for_comparison(result_a, method, n_act);
            u_b = embed_u_for_comparison(result_b, method, n_act);
            du = max(abs(u_a(:, common_ok) - u_b(:, common_ok)), [], 'all');
            dy = max(abs(data_a.y_achieved(:, common_ok) - data_b.y_achieved(:, common_ok)), [], 'all');
            inactive_a = inactive_u_max_abs(result_a, method, n_act, common_ok);
            inactive_b = inactive_u_max_abs(result_b, method, n_act, common_ok);
            [~, worst_linear_idx] = max(abs(u_a(:, common_ok) - u_b(:, common_ok)), [], 'all', 'linear');
            common_idx = find(common_ok);
            [~, worst_common_col] = ind2sub(size(u_a(:, common_ok)), worst_linear_idx);
            worst_sample = common_idx(worst_common_col);

            fprintf('[%-8s] samples=%d, max|du|=%.6g, max|d(Bu)|=%.6g, max inactive A=%.6g, max inactive B=%.6g, worst t=%.6g s\n', ...
                    method, nnz(common_ok), du, dy, inactive_a, inactive_b, t(worst_sample));
        end
    end
end

function idx = find_result_index_by_name(results, name)
    idx = -1;
    for i = 1:numel(results)
        if strcmp(results{i}.name, name)
            idx = i;
            return;
        end
    end
end

function n = get_plot_num_outputs(result)
    if isfield(result, 'plot_num_outputs') && ~isempty(result.plot_num_outputs)
        n = result.plot_num_outputs;
    else
        n = size(result.B, 2);
    end
end

function u_embed = embed_u_for_comparison(result, method, n_act)
    data = result.alloc.(method);
    u_embed = zeros(n_act, size(data.u, 2));

    if size(data.u, 1) == n_act
        u_embed = data.u;
        return;
    end

    if isfield(result, 'plot_active_idx') && numel(result.plot_active_idx) == size(data.u, 1)
        active_idx = result.plot_active_idx;
    else
        active_idx = 1:size(data.u, 1);
    end

    u_embed(active_idx, :) = data.u;
end

function value = inactive_u_max_abs(result, method, n_act, sample_mask)
    active = false(n_act, 1);

    if isfield(result, 'plot_active_idx') && ~isempty(result.plot_active_idx)
        active(result.plot_active_idx) = true;
    else
        active(1:size(result.alloc.(method).u, 1)) = true;
    end

    inactive = ~active;

    if ~any(inactive)
        value = 0;
        return;
    end

    u_embed = embed_u_for_comparison(result, method, n_act);
    value = max(abs(u_embed(inactive, sample_mask)), [], 'all');
end

function plot_reachable_sets_for_results(results)
    add_qcat_path_if_present();

    for case_idx = 1:numel(results)
        result = results{case_idx};

        [B_view, umin_view, umax_view, active_idx] = active_reachable_view_model(result.B, result.umin, result.umax);
        plim = [umin_view umax_view];

        figure('Name', ['Reachable set - ' result.name]);

        if exist('vview', 'file') ~= 2
            error('test:MissingQCAT', 'QCAT vview was not found on the MATLAB path.');
        end

        ratio = vview(B_view, plim, pinv(B_view));
        title(sprintf('%s reachable set, active u=[%s], ratio = %.6g', ...
              result.name, strtrim(sprintf('%d ', active_idx)), ratio), 'Interpreter', 'none');
    end
end

function [B_view, umin_view, umax_view, active_idx] = active_reachable_view_model(B, umin, umax)
    limit_active = abs(umax - umin) > 1e-9;
    column_active = vecnorm(B, 2, 1)' > 1e-12;
    active_idx = find(limit_active & column_active);

    if isempty(active_idx)
        active_idx = 1:size(B, 2);
    end

    B_view = B(:, active_idx);
    umin_view = umin(active_idx);
    umax_view = umax(active_idx);
end

function add_qcat_path_if_present()
    qcat_paths = {'QCAT/qcat', 'control_allocation_lib/qcat/QCAT/qcat'};

    for i = 1:numel(qcat_paths)
        if isfolder(qcat_paths{i})
            addpath(genpath(qcat_paths{i}));
        end
    end
end
