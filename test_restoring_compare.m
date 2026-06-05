clear;
clc;
addpath(genpath(pwd));

%% Compare restoring_cpp against the original restoring.m
% restoring_cpp is the C++-aligned version used by test.m.
% restoring is the original MATLAB version based on pinv([B;u']).

load 'input.mat';
if isfile('input.csv')
    v = readmatrix('input.csv')';
end

dt = mean(controls_delta_t_s);
t_full = 0:dt:dt*(len_command_px4-1);
test_time_window_s = [20 30];
sample_indices = find(t_full >= test_time_window_s(1) & t_full <= test_time_window_s(2));
v_test = v(:, sample_indices);
t = t_full(sample_indices);

cases = {make_aircraft_4(), make_aircraft_6()};
restorers = { ...
    struct('name', 'restoring_cpp', 'fn', @restoring_cpp), ...
    struct('name', 'restoring', 'fn', @restoring), ...
    struct('name', 'restoring_opt', 'fn', @restoring_opt), ...
    struct('name', 'restoring_opt1', 'fn', @restoring_opt1), ...
    struct('name', 'restoring_opt2', 'fn', @restoring_opt2), ...
    struct('name', 'restoring_opt3', 'fn', @restoring_opt3), ...
    struct('name', 'restoring_return_k', 'fn', @restoring_return_k), ...
    struct('name', 'restoring_robustlyB__return_k', 'fn', @restoring_robustlyB__return_k) ...
};

for case_idx = 1:numel(cases)
    aircraft = cases{case_idx};
    B = aircraft.B;
    umin = aircraft.umin;
    umax = aircraft.umax;
    [k, m] = size(B);
    N = size(v_test, 2);
    m_higher = zeros(k, 1);

    raw = zeros(m, N);
    outputs = struct();
    failures = struct();
    for ridx = 1:numel(restorers)
        outputs.(restorers{ridx}.name) = NaN(m, N);
        failures.(restorers{ridx}.name) = 0;
    end

    for idx = 1:N
        [u, ~, ~] = DP_LPCA_copy(m_higher, v_test(:, idx), B, umin, umax, 100);
        u = min(max(u, umin), umax);
        raw(:, idx) = u;
        for ridx = 1:numel(restorers)
            name = restorers{ridx}.name;
            try
                outputs.(name)(:, idx) = restorers{ridx}.fn(B, u, umin, umax);
            catch ME
                failures.(name) = failures.(name) + 1;
                if failures.(name) == 1
                    fprintf('[%s %s] first failure at sample %d, time %.6g: %s\n', ...
                        aircraft.name, name, idx, t(idx), ME.message);
                end
            end
        end
    end

    ref = outputs.restoring_cpp;
    fprintf('\n[%s restoring methods vs restoring_cpp, %.6g-%.6g s]\n', aircraft.name, t(1), t(end));
    fprintf('%-32s %14s %14s %14s %10s\n', 'restorer', 'max |du|', 'max |Bdu|', 'max |u-raw|', 'failures');
    for ridx = 1:numel(restorers)
        name = restorers{ridx}.name;
        x = outputs.(name);
        valid = all(isfinite(x), 1);
        if any(valid)
            du = max(abs(x(:, valid) - ref(:, valid)), [], 'all');
            dBu = max(abs(B * (x(:, valid) - ref(:, valid))), [], 'all');
            move = max(abs(x(:, valid) - raw(:, valid)), [], 'all');
        else
            du = NaN;
            dBu = NaN;
            move = NaN;
        end
        fprintf('%-32s %14.6g %14.6g %14.6g %10d\n', name, du, dBu, move, failures.(name));
    end
end

function aircraft = make_aircraft_4()
    l1 = 0.167;
    l2 = 0.069;
    k_v = 3;
    I = diag([0.01149, 0.01153, 0.00487]);
    B = I \ [-l1 0 l1 0;
             0 -l1 0 l1;
             l2 l2 l2 l2] * k_v;
    ulim = 0.3491;
    aircraft.name = '4-effector';
    aircraft.B = B;
    aircraft.umin = ones(4, 1) * -ulim;
    aircraft.umax = ones(4, 1) * ulim;
end

function aircraft = make_aircraft_6()
    l1 = 0.2998;
    l2 = 0.0664;
    I = diag([0.05, 0.04, 0.013]);
    d = 60*pi/180;
    B = I \ [-l1, -l1*cos(d),  l1*cos(d),  l1,  l1*cos(d), -l1*cos(d);
              0,   l1*sin(d),  l1*sin(d),  0,  -l1*sin(d), -l1*sin(d);
              l2,  l2,         l2,         l2,  l2,         l2];
    ulim = 40*pi/180;
    aircraft.name = '6-effector';
    aircraft.B = B;
    aircraft.umin = ones(6, 1) * -ulim;
    aircraft.umax = ones(6, 1) * ulim;
end
