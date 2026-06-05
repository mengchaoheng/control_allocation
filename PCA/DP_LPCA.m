function [u, errout, lambda] = DP_LPCA(yd, B, uMin, uMax, itlim, varargin)
% Direction Preserving Control Allocation Linear Program.
%
% Top-level version of the book DP_LPCA routine.  The problem construction
% matches control_allocation_lib/aircraft-control-allocation-book-simulation,
% and the simplex backend is the original book/MATLAB simplxuprevsol.m.
% Extra varargin inputs are accepted only for legacy scripts and are ignored.

errout = 0;
lambda = 0;
[n, m] = size(B);

if all(abs(yd) < eps)
    errout = -1;
    u = zeros(m, 1);
    return;
end

A = [B -yd];
b = -B * uMin;
c = [zeros(m, 1); -1];
h = [uMax - uMin; 1];

sb = 2 * (b > 0) - 1;
Ai = [A diag(sb)];
ci = [zeros(m + 1, 1); ones(n, 1)];
inBi = m + 2:m + n + 1;
ei = true(m + n + 1, 1);
hi = [h; 2 * abs(b)];

[y1, inB1, e1, itlim, errsimp] = call_simplex(Ai, ci', b, inBi, hi, ei, n, m + n + 1, itlim);

if itlim <= 0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1 > (m + 1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~= 0
    xout = zeros(m + 1, 1);
    indv = inB1 <= (m + 1);
    xout(inB1(indv)) = y1(indv);
    xout(~e1(1:m + 1)) = -xout(~e1(1:m + 1)) + h(~e1(1:m + 1));
else
    [y2, inB2, e2, itlim, errsimp] = call_simplex(A, c', b, inB1, h, e1(1:m + 1), n, m + 1, itlim);

    xout = zeros(m + 1, 1);
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

u = xout(1:m) + uMin;
lambda = xout(m + 1);

end

function [y, inB, e, itlim, errsimp] = call_simplex(A, c, b, inB, h, e, m, n, itlim)
    [y, inB, e, itlim, errsimp] = simplxuprevsol(A, c, b, inB, h, e, m, n, itlim);
end
