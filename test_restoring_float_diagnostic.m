clear all;
close all;
clc;
addpath(genpath(pwd));

if ~isfile('test_results.mat')
    error('test_restoring_float_diagnostic:MissingResults', ...
          'Run test.m first so test_results.mat contains cpp_dir raw outputs.');
end

data = load('test_results.mat');
result = find_result_by_name(data.results, '15006_SHC09');
if isempty(result)
    error('test_restoring_float_diagnostic:MissingCase', ...
          'test_results.mat does not contain case 15006_SHC09.');
end
if ~isfield(result, 'alloc') || ~isfield(result.alloc, 'cpp_dir')
    error('test_restoring_float_diagnostic:MissingCppDir', ...
          'test_results.mat does not contain result.alloc.cpp_dir.');
end

B = result.B;
umin = result.umin;
umax = result.umax;
u_raw_all = result.alloc.cpp_dir.u_raw;
u_cpp_all = result.alloc.cpp_dir.u;
u_matlab_all = result.alloc.dir.u;
n = size(u_raw_all, 2);

u_double = zeros(size(u_raw_all));
u_float_normal = zeros(size(u_raw_all));
u_float_mgs = zeros(size(u_raw_all));
u_float_normal_relaxed = zeros(size(u_raw_all));

diag_normal = repmat(empty_diag(), 1, n);
diag_mgs = repmat(empty_diag(), 1, n);
diag_normal_relaxed = repmat(empty_diag(), 1, n);

for idx = 1:n
    u_raw = u_raw_all(:, idx);
    u_double(:, idx) = restoring_cpp(B, u_raw, umin, umax);
    [u_float_normal(:, idx), diag_normal(idx)] = restoring_float_normal_equation(B, u_raw, umin, umax);
    [u_float_normal_relaxed(:, idx), diag_normal_relaxed(idx)] = restoring_float_normal_equation(B, u_raw, umin, umax, 1e-4);
    [u_float_mgs(:, idx), diag_mgs(idx)] = restoring_float_mgs_projection(B, u_raw, umin, umax);
end

fprintf('\n=== restoring float diagnostic: 15006_SHC09 cpp_dir raw input ===\n');
fprintf('samples = %d\n', n);
fprintf('current C++ restored vs MATLAB restoring_cpp(cpp raw): %.12g\n', max(abs(u_cpp_all - u_double), [], 'all'));
fprintf('float normal-equation restoring vs MATLAB restoring_cpp(cpp raw): %.12g\n', max(abs(u_float_normal - u_double), [], 'all'));
fprintf('float normal-equation restoring, residual tol=1e-4, vs MATLAB restoring_cpp(cpp raw): %.12g\n', max(abs(u_float_normal_relaxed - u_double), [], 'all'));
fprintf('float MGS restoring vs MATLAB restoring_cpp(cpp raw): %.12g\n', max(abs(u_float_mgs - u_double), [], 'all'));
fprintf('MATLAB dir restored vs MATLAB restoring_cpp(cpp raw): %.12g\n', max(abs(u_matlab_all - u_double), [], 'all'));

[~, worst_normal_linear] = max(abs(u_float_normal - u_double), [], 'all', 'linear');
[~, worst_mgs_linear] = max(abs(u_float_mgs - u_double), [], 'all', 'linear');
[worst_normal_u, worst_normal_sample] = ind2sub(size(u_float_normal), worst_normal_linear);
[worst_mgs_u, worst_mgs_sample] = ind2sub(size(u_float_mgs), worst_mgs_linear);

print_sample('normal-equation worst', worst_normal_sample, worst_normal_u, data, result, ...
             u_raw_all, u_double, u_float_normal, diag_normal);
print_sample('normal-equation relaxed same sample', worst_normal_sample, worst_normal_u, data, result, ...
             u_raw_all, u_double, u_float_normal_relaxed, diag_normal_relaxed);
print_sample('MGS worst', worst_mgs_sample, worst_mgs_u, data, result, ...
             u_raw_all, u_double, u_float_mgs, diag_mgs);

function d = empty_diag()
    d.null_norm2 = NaN;
    d.residual_norm = NaN;
    d.k_opt = NaN;
    d.k_max = NaN;
    d.k_used = NaN;
    d.return_reason = "";
end

function [u_rest, d] = restoring_float_normal_equation(B, u, uMin, uMax, residual_tol_value)
    tol = single(1e-5);
    if nargin < 5
        residual_tol = tol;
    else
        residual_tol = single(residual_tol_value);
    end
    a = single(-2);
    d = empty_diag();

    Bf = single(B);
    uf = single(u);
    uMinf = single(uMin);
    uMaxf = single(uMax);

    if norm(uf) < tol
        u_rest = double(uf);
        d.return_reason = "zero-u";
        return;
    end

    achieved = Bf * uf;
    row_solution = (Bf * Bf') \ achieved;
    u_pseudo = Bf' * row_solution;
    null_component = uf - u_pseudo;
    null_norm2 = null_component' * null_component;
    d.null_norm2 = double(null_norm2);

    if null_norm2 < tol^2
        u_rest = double(uf);
        d.return_reason = "zero-null";
        return;
    end

    u_null = a * null_component / null_norm2;
    residual_norm = norm(Bf * u_null);
    d.residual_norm = double(residual_norm);
    if residual_norm > residual_tol
        u_rest = double(uf);
        d.return_reason = "residual";
        return;
    end

    u_null_norm2 = u_null' * u_null;
    if u_null_norm2 < tol^2
        u_rest = double(uf);
        d.return_reason = "zero-unull";
        return;
    end

    k_opt = -a / u_null_norm2;
    [k_max, k_used] = restoring_step_limit(uf, u_null, uMinf, uMaxf, k_opt, tol);
    u_rest = double(uf + k_used * u_null);
    d.k_opt = double(k_opt);
    d.k_max = double(k_max);
    d.k_used = double(k_used);
    d.return_reason = "restored";
end

function [u_rest, d] = restoring_float_mgs_projection(B, u, uMin, uMax)
    tol = single(1e-5);
    a = single(-2);
    d = empty_diag();

    Bf = single(B);
    uf = single(u);
    uMinf = single(uMin);
    uMaxf = single(uMax);
    [k, m] = size(Bf);

    if norm(uf) < tol
        u_rest = double(uf);
        d.return_reason = "zero-u";
        return;
    end

    Q = zeros(k, m, 'single');
    q_count = 0;
    for row = 1:k
        v = Bf(row, :);
        for qi = 1:q_count
            v = v - dot(v, Q(qi, :)) * Q(qi, :);
        end
        v_norm = norm(v);
        if v_norm > tol
            q_count = q_count + 1;
            Q(q_count, :) = v / v_norm;
        end
    end

    u_pseudo = zeros(m, 1, 'single');
    for qi = 1:q_count
        q = Q(qi, :)';
        u_pseudo = u_pseudo + q * (q' * uf);
    end

    null_component = uf - u_pseudo;
    null_norm2 = null_component' * null_component;
    d.null_norm2 = double(null_norm2);
    if null_norm2 < tol^2
        u_rest = double(uf);
        d.return_reason = "zero-null";
        return;
    end

    u_null = a * null_component / null_norm2;
    residual_norm = norm(Bf * u_null);
    d.residual_norm = double(residual_norm);
    if residual_norm > tol
        u_rest = double(uf);
        d.return_reason = "residual";
        return;
    end

    u_null_norm2 = u_null' * u_null;
    if u_null_norm2 < tol^2
        u_rest = double(uf);
        d.return_reason = "zero-unull";
        return;
    end

    k_opt = -a / u_null_norm2;
    [k_max, k_used] = restoring_step_limit(uf, u_null, uMinf, uMaxf, k_opt, tol);
    u_rest = double(uf + k_used * u_null);
    d.k_opt = double(k_opt);
    d.k_max = double(k_max);
    d.k_used = double(k_used);
    d.return_reason = "restored";
end

function [k_max, k_used] = restoring_step_limit(u, u_null, uMin, uMax, k_opt, tol)
    k_max = single(Inf);
    for i = 1:numel(u)
        if abs(u_null(i)) < tol
            continue;
        end
        if u_null(i) > 0
            candidate = (uMax(i) - u(i)) / u_null(i);
        else
            candidate = (uMin(i) - u(i)) / u_null(i);
        end
        if candidate < k_max
            k_max = candidate;
        end
    end
    k_used = min(k_max, k_opt);
end

function print_sample(label, sample_idx, u_idx, data, result, u_raw_all, u_double, u_test, diag_data)
    fprintf('\n[%s]\n', label);
    if isfield(data, 't') && sample_idx <= numel(data.t)
        fprintf('sample=%d, t=%.9g s, worst_u=%d\n', sample_idx, data.t(sample_idx), u_idx);
    else
        fprintf('sample=%d, worst_u=%d\n', sample_idx, u_idx);
    end
    fprintf('u_raw      = ['); fprintf(' %.12g', u_raw_all(:, sample_idx)); fprintf(' ]^T\n');
    fprintf('u_double   = ['); fprintf(' %.12g', u_double(:, sample_idx)); fprintf(' ]^T\n');
    fprintf('u_test     = ['); fprintf(' %.12g', u_test(:, sample_idx)); fprintf(' ]^T\n');
    fprintf('max |du|   = %.12g\n', max(abs(u_test(:, sample_idx) - u_double(:, sample_idx))));
    fprintf('B*u_double = ['); fprintf(' %.12g', result.B * u_double(:, sample_idx)); fprintf(' ]^T\n');
    fprintf('B*u_test   = ['); fprintf(' %.12g', result.B * u_test(:, sample_idx)); fprintf(' ]^T\n');
    d = diag_data(sample_idx);
    fprintf('null_norm2=%.12g, residual_norm=%.12g, k_opt=%.12g, k_max=%.12g, k_used=%.12g, reason=%s\n', ...
            d.null_norm2, d.residual_norm, d.k_opt, d.k_max, d.k_used, d.return_reason);
end

function result = find_result_by_name(results, name)
    result = [];
    for idx = 1:numel(results)
        if strcmp(results{idx}.name, name)
            result = results{idx};
            return;
        end
    end
end
