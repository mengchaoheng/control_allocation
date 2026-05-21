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
    % // B=K*P, K=I\diag([l1 l1 l2])k, P=[-1     0     1     0; 0    -1     0     1; 1     1     1     1];   K=diag([ k*l1/I_x  k*l1/I_y  k*l2/I_z  ])
    %     // B^{\dagger} = P^\top K K^{-1} (P P^\top)^{-1} K^{-1} = P^\top (P P^\top)^{-1} K^{-1}=P^{\dagger} K^{-1}
    % // P^{\dagger}=[-0.5000   -0.0000    0.2500;0   -0.5000    0.2500;0.5000   -0.0000    0.2500;0    0.5000    0.2500]
    % // B^{\dagger}=P^{\dagger} K^{-1}=[-0.5000   -0.0000    0.2500;0   -0.5000    0.2500;0.5000   -0.0000    0.2500;0    0.5000    0.2500]*diag([ I_x/(k*l1)  I_y/(k*l1)  I_z/(k*l2)  ])

load 'input.mat'; % provides v, len_command_px4, controls_delta_t_s, ...
if isfile('input.csv')
    v = readmatrix('input.csv')';
end

[~, N] = size(v);
dt = mean(controls_delta_t_s);
t_full = 0:dt:dt*(len_command_px4-1);

% Use [] for the full input, or a time window such as [20 30] to speed up
% focused comparisons.  The C++ CSV outputs are sliced with the same
% original sample indices, so reports still compare the same flight segment.
test_time_window_s = [20 30];
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
% To test another B, append a case here:
% aircraft_cases{end+1} = make_aircraft_from_matrix('my-B', B, umin, umax, '');
% If cpp_tag is '4' or '6', test.m loads results/cpp_outputs/output_cpp_<tag>_*.csv for C++ comparison.
% If cpp_tag is '', only MATLAB outputs are simulated and saved/plotted.

%% Tie-break settings shared with alloc_cpp/src/ControlAllocation/ControlAllocation.h
tie_opts.tie_rel_tol = 1e-5;
tie_opts.tie_abs_tol = 1e-6;
tie_opts.zero_tie_abs_tol = 3e-5;

fprintf('MATLAB path: PCA/DP_LPCA, PCA/DPscaled_LPCA, PCA/DP_LPCA_prio + restoring_cpp\n');
fprintf('C++ PCA path: DP_LPCA / DPscaled_LPCA / DP_LPCA_prio + restoring() aligned to restoring_cpp\n');

%% Simulate flight process
tic;
results = cell(size(aircraft_cases));
for case_idx = 1:numel(aircraft_cases)
    results{case_idx} = simulate_flight_process(aircraft_cases{case_idx}, v_test, tie_opts);
end
elapsed_time = toc;
fprintf('MATLAB reference execution time: %.2f 秒\n', elapsed_time);
fprintf('Test sample window: %.6g to %.6g s (%d samples)\n', t(1), t(end), len_command_px4);

%% C++ outputs
ensure_cpp_outputs();

for case_idx = 1:numel(results)
    results{case_idx} = load_cpp_outputs_for_case(results{case_idx}, test_sample_indices);
end

%% Reports
for case_idx = 1:numel(results)
    report_case(results{case_idx}, command_px4, len_command_px4, t);
end

%% Save once, then rerun only plot_test_results.m when changing plots
result4 = results{1};
result6 = results{2};
save('test_results.mat', 'results', 'result4', 'result6', 'command_px4', 'len_command_px4', 't', 'tie_opts', 'test_time_window_s', 'test_sample_indices');
plot_test_results('test_results.mat');

function aircraft = make_aircraft_4()
    % Four-effector ducted-fan model.  This form keeps the geometry and
    % inertia parameters visible; it is numerically equal to alloc_cpp/test/main.cpp.
    %
    % 行满秩时 pinv(B) = The Moore-Penrose Pseudo-inverse
    % If B is full row rank, The Moore-Penrose Pseudo-inverse
    % B^+ = B^T (B B^T)^{-1}, since
    % B = K*P,
    % K = I\diag([l1 l1 l2])*k,
    % P = [-1     0     1     0;
    %       0    -1     0     1;
    %       1     1     1     1];
    % K = diag([ k*l1/I_x  k*l1/I_y  k*l2/I_z ])
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
    %            = P^{dagger} * diag([ I_x/(k*l1)  I_y/(k*l1)  I_z/(k*l2) ])

    l1 = 0.167;
    l2 = 0.069;
    k_v = 3;
    I_x = 0.01149;
    I_y = 0.01153;
    I_z = 0.00487;
    I = diag([I_x, I_y, I_z]);

    %=============================4==================================
    B = I \ [-l1     0       l1     0;
              0      -l1     0      l1;
              l2     l2      l2     l2] * k_v;
    % B = I\diag([l1 l1 l2])[-1     0     1     0;
    %                         0    -1     0     1;
    %                         1     1     1     1]*k_v;
    %
    % B^+ = B^T (B B^T)^{-1}
    %    1X
    % 4     2Y
    %    3

    ulim = 0.3491;
    aircraft.name = '4-effector';
    aircraft.cpp_tag = '4';
    aircraft.B = B;
    aircraft.umin = ones(4, 1) * -ulim;
    aircraft.umax = ones(4, 1) * ulim;
end

function aircraft = make_aircraft_6()
    % Six-effector model, same physical construction as alloc_cpp/test/main.cpp.
    l1 = 0.2998;
    l2 = 0.0664;
    k_omega2force = 1.0;
    I_x = 0.05;
    I_y = 0.04;
    I_z = 0.013;
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

    ulim = 40*pi/180;
    aircraft.name = '6-effector';
    aircraft.cpp_tag = '6';
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

function result = simulate_flight_process(aircraft, v, tie_opts)
    B = aircraft.B;
    umin = aircraft.umin;
    umax = aircraft.umax;
    [k, m] = size(B);
    [~, N] = size(v);
    m_higher = zeros(k, 1);

    result = aircraft;
    result.matlab_dp_raw = zeros(m, N);
    result.matlab_dpscaled_raw = zeros(m, N);
    result.matlab_prio_raw = zeros(m, N);
    result.matlab_dp = zeros(m, N);
    result.matlab_dpscaled = zeros(m, N);
    result.matlab_prio = zeros(m, N);

    %% setup function of allocation lib
    % ========
    %% setup ACA
    global NumU
    NumU=m;
    LPmethod=2; % LPmethod should be an integer between 0 and 5. when LPmethod=2 set upper of lambda to Inf can't save this method!!! but big number is the same as that method based linprog
    % DPscaled_LPCA的结果和restoring接近但是有细微区别，对lambda限制在0-1之间，有助于优先级的理论推导。
    INDX=ones(1,m);  % active effectors
    IN_MAT = [B     zeros(k,1)
              umin' 0
              umax' 0
              INDX  LPmethod];
    % for frame-wise
    % u_0=ones(m,1)*20*pi/180;
    % u_0=[0.0122;
    %      0.0122;
    %      0.3491;
    %      0.3491];
    % u_0=[0.0;
    %      0.0;
    %      0.0;
    %      0.0];
    % IN_MAT1 = [B     zeros(k,1)
    %           (umin-u_0)' 0
    %           (umax-u_0)' 0
    %           INDX  LPmethod];
    %% setup qcat. just wls_alloc and not a Hotstart setting here, use test_qcat.m for more test,
    Wv   = eye(k);     % QP allocation
    Wu   = eye(m);
    ud   = zeros(m,1);
    gam  = 1e6;        % weight
    u0 = (umin+umax)/2;
    W0 = zeros(m,1);
    imax = 100;        % no of iterations
    %
    % ========
    %%
    u=zeros(m,1);
    x_LPwrap=zeros(m,N);
    x_LPwrap1=zeros(m,N);
    x_PCA=zeros(m,N);
    x_LPwrap_incre=zeros(m,N);
    x_allocator_dir_LPwrap_4=zeros(m,N);
    x_CGIwrapp=zeros(m,N);
    x_DAwrap=zeros(m,N);
    x_VJAwrap=zeros(m,N);
    x_inv=zeros(m,N);
    x_wls=zeros(m,N);
    x_wls_gen=zeros(m,N);
    x_dir_alloc_linprog=zeros(m,N);
    x_dir_alloc_linprog_re=zeros(m,N);
    x_dir_alloc_linprog_re_bound=zeros(m,N);
    x_use_LP_lib=zeros(m,N);

    %% simulate flight process
    for idx = 1:N  % or x:N for debug
        command = v(:, idx);
        IN_MAT(1:k,end) = command + m_higher; %[ 36.8125; 0;92.9776];%

        [u, ~, ~] = DP_LPCA(command, B, umin, umax, 100, tie_opts);
        u = min(max(u, umin), umax);
        result.matlab_dp_raw(:, idx) = u;
        result.matlab_dp(:, idx) = restoring_cpp(B, u, umin, umax);

        [u, ~, ~, ~] = DPscaled_LPCA(command, B, umin, umax, 100, tie_opts);
        u = min(max(u, umin), umax);
        result.matlab_dpscaled_raw(:, idx) = u;
        result.matlab_dpscaled(:, idx) = restoring_cpp(B, u, umin, umax);

        [u, ~, ~] = DP_LPCA_prio(m_higher, command, B, umin, umax, 100, tie_opts);
        u = min(max(u, umin), umax);
        result.matlab_prio_raw(:, idx) = u;
        result.matlab_prio(:, idx) = restoring_cpp(B, u, umin, umax);

        % Add extra allocators here when needed. Keep each output as an
        % m-by-N field in result, e.g. result.wls(:, idx), so
        % plot_test_results.m can redraw it after one run.

        % IN_MAT(end,end) = 2; % LPmethod=2: original DP_LPCA in LPwrap.m
        % u = LPwrap(IN_MAT); % function of ACA lib
        % u = min(max(u, umin), umax);
        % result.LPwrap_DP_LPCA_raw(:,idx) = u;
        % result.LPwrap_DP_LPCA(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % IN_MAT(end,end) = 3; % LPmethod=3: original DPscaled_LPCA in LPwrap.m
        % u = LPwrap(IN_MAT); % function of ACA lib
        % u = min(max(u, umin), umax);
        % result.LPwrap_DPscaled_LPCA_raw(:,idx) = u;
        % result.LPwrap_DPscaled_LPCA(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % IN_MAT(end,end) = LPmethod;
        % u = LPwrap(IN_MAT); % function of ACA lib
        % u = min(max(u, umin), umax);
        % result.LPwrap_raw(:,idx) = u;
        % result.LPwrap(:,idx) = restoring_cpp(B, u, umin, umax);

        % [u, ~, ~] = DP_LPCA_prio(m_higher, command, B, umin, umax, 100, tie_opts);
        % u = min(max(u, umin), umax);
        % result.PCA_raw(:,idx) = u;
        % result.PCA(:,idx) = restoring_cpp(B, u, umin, umax);

        % u = LPwrap(IN_MAT1); % incremental form.
        % The control constraint δ ≤ δ ≤ δ must contain the origin, i.e., δ = 0
        % must be a feasible control input.
        % In order word, 0 have to be a feasible solution. When add vel contraint to s.t.
        % and if optimization variables is delta_u, then 0 is a feasible
        % solution, else if optimization variables is u, then 0 is not a
        % feasible solution, so the LP Not working.
        % u = min(max(u, umin-u_0), umax-u_0) + u_0;
        % result.LPwrap_incre_raw(:,idx) = u;
        % result.LPwrap_incre(:,idx) = restoring_cpp(B, u, umin, umax);

        % u = CGIwrap(IN_MAT);
        % u = min(max(u, umin), umax);
        % result.CGIwrap_raw(:,idx) = u;
        % result.CGIwrap(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % u = DAwrap(IN_MAT);
        % u = min(max(u, umin), umax);
        % result.DAwrap_raw(:,idx) = u;
        % result.DAwrap(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % u = VJAwrap(IN_MAT);
        % u = min(max(u, umin), umax);
        % result.VJAwrap_raw(:,idx) = u;
        % result.VJAwrap(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % u = pinv(B) * command;
        % u = min(max(u, umin), umax);
        % result.inv_raw(:,idx) = u;
        % result.inv(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % [u, ~, ~] = wls_alloc(B, command, umin, umax, Wv, Wu, ud, gam, u0, W0, imax);
        % u = min(max(u, umin), umax);
        % result.wls_raw(:,idx) = u;
        % result.wls(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % u = wls_alloc_gen(B, command, umin, umax, eye(k), eye(m), zeros(m,1), 1e6, zeros(m,1), zeros(m,1), 100, 4);
        % u = min(max(u, umin), umax);
        % result.wls_gen_raw(:,idx) = u;
        % result.wls_gen(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % [u, ~] = dir_alloc_linprog(B, command, umin, umax, 1e4); % LPmethod=2 and lam=1 of dir_alloc_linprog is lager but similar
        % u = min(max(u, umin), umax);
        % result.dir_alloc_linprog_raw(:,idx) = u;
        % result.dir_alloc_linprog(:,idx) = restoring_cpp(B, u, umin, umax);
        %
        % [u, ~] = dir_alloc_linprog_re(B, command, umin, umax);
        % u = min(max(u, umin), umax);
        % result.dir_alloc_linprog_re_raw(:,idx) = u;
        % result.dir_alloc_linprog_re(:,idx) = restoring_cpp(B, u, umin, umax);

        % [u, ~] = dir_alloc_linprog_re_bound(B, command, umin, umax, 1e4);% the
        % same as dir_alloc_linprog for any lam >=1, lam have to be >1 when use
        % linprog, that will be the same as LPmethod=3
        % u = min(max(u, umin), umax);
        % result.dir_alloc_linprog_re_bound_raw(:,idx) = u;
        % result.dir_alloc_linprog_re_bound(:,idx) = restoring_cpp(B, u, umin, umax);

        % [u, ~] = use_LP_lib(B, command, umin, umax); % ToDo: use the LP lib
        % u = min(max(u, umin), umax);
        % result.use_LP_lib_raw(:,idx) = u;
        % result.use_LP_lib(:,idx) = restoring_cpp(B, u, umin, umax);

        % [u, ~, ~] = allocator_dir_LPwrap_4(single(B), single(command), single(umin), single(umax)); % ToDo: 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.914283e-09。
        % u = min(max(u, umin), umax);
        % result.allocator_dir_LPwrap_4_raw(:,idx) = u;
        % result.allocator_dir_LPwrap_4(:,idx) = restoring_cpp(B, u, umin, umax);
    end
end

function ensure_cpp_outputs()
    cpp_output_dir = fullfile('results', 'cpp_outputs');
    cpp_output_files = [fullfile(cpp_output_dir, "output_cpp_4_dp.csv"), ...
                        fullfile(cpp_output_dir, "output_cpp_4_dpscaled.csv"), ...
                        fullfile(cpp_output_dir, "output_cpp_4_prio.csv"), ...
                        fullfile(cpp_output_dir, "output_cpp_6_dp.csv"), ...
                        fullfile(cpp_output_dir, "output_cpp_6_dpscaled.csv"), ...
                        fullfile(cpp_output_dir, "output_cpp_6_prio.csv")];
    if all(isfile(cpp_output_files))
        return;
    end
    if ~isfile('./alloc_cpp/build/main')
        error('test:MissingCppBinary', 'Build alloc_cpp/build/main before running test.m.');
    end
    status = system('./alloc_cpp/build/main');
    if status ~= 0
        error('test:CppRunFailed', 'Failed to run alloc_cpp/build/main.');
    end
end

function result = load_cpp_outputs_for_case(result, sample_indices)
    if ~isfield(result, 'cpp_tag') || isempty(result.cpp_tag)
        return;
    end

    tag = result.cpp_tag;
    cpp_output_dir = fullfile('results', 'cpp_outputs');
    dp_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_dp.csv', tag));
    dpscaled_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_dpscaled.csv', tag));
    prio_file = fullfile(cpp_output_dir, sprintf('output_cpp_%s_prio.csv', tag));

    if isfile(dp_file)
        result.cpp_dp = slice_cpp_output(readmatrix(dp_file)', sample_indices);
    end
    if isfile(dpscaled_file)
        result.cpp_dpscaled = slice_cpp_output(readmatrix(dpscaled_file)', sample_indices);
    end
    if isfile(prio_file)
        result.cpp_prio = slice_cpp_output(readmatrix(prio_file)', sample_indices);
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

function report_case(result, command_px4, len_command_px4, t)
    if isfield(result, 'cpp_dp')
        report_cpp_matlab_case([result.name ' DP_LPCA'], result.B, result.cpp_dp, result.matlab_dp, command_px4, len_command_px4, t);
    else
        report_matlab_case([result.name ' DP_LPCA'], result.B, result.matlab_dp, command_px4, len_command_px4);
    end

    if isfield(result, 'cpp_dpscaled')
        report_cpp_matlab_case([result.name ' DPscaled_LPCA'], result.B, result.cpp_dpscaled, result.matlab_dpscaled, command_px4, len_command_px4, t);
    else
        report_matlab_case([result.name ' DPscaled_LPCA'], result.B, result.matlab_dpscaled, command_px4, len_command_px4);
    end

    if isfield(result, 'cpp_prio')
        report_cpp_matlab_case([result.name ' DP_LPCA_prio'], result.B, result.cpp_prio, result.matlab_prio, command_px4, len_command_px4, t);
    else
        report_matlab_case([result.name ' DP_LPCA_prio'], result.B, result.matlab_prio, command_px4, len_command_px4);
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
