function [u, itlim, errout, rho] = DPscaled_LPCA(yd, B, uMin, uMax, itlim, varargin)
% Direction Preserving Control Allocation Linear Program, reduced form.
%
% Top-level version of the book DPscaled_LPCA routine.  The problem
% construction matches control_allocation_lib/aircraft-control-allocation-
% book-simulation, and the simplex backend is the original book/MATLAB
% simplxuprevsol.m.
% Extra varargin inputs are accepted only for legacy scripts and are ignored.

errout = 0;
rho = 0;
[n, m] = size(B);

[my, iy] = max(abs(yd));
if my < eps
    errout = -1;
    u = zeros(m, 1);
    return;
end

Bt = B([iy setdiff(1:n, iy)], :);
ydt = yd([iy setdiff(1:n, iy)]);
ydt(2:3) = ydt([3 2]);
Bt([2 3], :) = Bt([3 2], :);

M = [ydt(2:n) -ydt(1) * eye(n - 1)];
A = M * Bt;
b = -A * uMin;
c = -Bt' * ydt;
h = uMax - uMin;

sb = 2 * (b > 0) - 1;
Ai = [A diag(sb)];
ci = [zeros(m, 1); ones(n - 1, 1)];
inBi = m + 1:m + n - 1;
ei = true(m + n - 1, 1);
hi = [h; 2 * abs(b)];

[y1, inB1, e1, itlim, errsimp] = call_simplex(Ai, ci', b, inBi, hi, ei, n - 1, m + n - 1, itlim);

if itlim <= 0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1 > m)
    errout = -2;
    disp('No Initial Feasible Solution found');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~= 0
    xout = zeros(m, 1);
    indv = inB1 <= m;
    xout(inB1(indv)) = y1(indv);
    xout(~e1(1:m)) = -xout(~e1(1:m)) + h(~e1(1:m));
else
    [y2, inB2, e2, itlim, errsimp] = call_simplex(A, c', b, inB1, h, e1(1:m), n - 1, m, itlim);

    xout = zeros(m, 1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2) + h(~e2);

    if itlim <= 0
        errout = 3;
        disp('Too Many Iterations Finding Final Solution');
    end
    if errsimp
        errout = 1;
        disp('Solver error');
    end
end

u = xout + uMin;
rho = ydt' * Bt * u / (ydt' * ydt);
if rho > 1
    u = u / rho;
end

end

function [y, inB, e, itlim, errsimp] = call_simplex(A, c, b, inB, h, e, m, n, itlim)
    [y, inB, e, itlim, errsimp] = simplxuprevsol(A, c, b, inB, h, e, m, n, itlim);
end
