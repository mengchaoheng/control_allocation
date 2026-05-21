function [u, errout, lambda] = DP_LPCA_copy(m_higher, yd, B, uMin, uMax, itlim, opts)
% Direction Preserving Control Allocation Linear Program.
% only change b = m_higher - B*uMin;
%
% This is the maintained MATLAB DP_LPCA_copy implementation.  It keeps the
% original DP_LPCA_copy interface and problem construction, but uses
% simplxuprevsol_tiebreak as the bounded revised simplex backend so MATLAB
% follows the same deterministic pivot/tie-break rules as the C++ allocator.
%
% function [u, errout, lambda] = DP_LPCA_copy(m_higher, yd, B, uMin, uMax, itlim, opts)
%
% Inputs:
%         m_higher [n] = higher-priority objective offset
%         yd [n]       = lower-priority desired objective direction
%         B [n,m]      = Control Effectiveness matrix
%         uMin[m,1]    = Lower bound for controls
%         uMax[m,1]    = Upper bound for controls
%         itlim        = Number of allowed iterations limit
%         opts         = optional tie-break settings for simplxuprevsol_tiebreak
%
% Outputs:
%         u[m,1]       = Control Solution
%         errout       = Error Status code
%         lambda       = direction-preserving scale factor
%
% Calls:
%         simplxuprevsol_tiebreak = Bounded Revised Simplex solver with
%         deterministic tie-breaks.

if nargin < 7
    opts = struct();
end

errout = 0;
lambda = 0;
[n, m] = size(B);

% Check to see if yd == 0.
if all(abs(yd) < eps)
    errout = -1;
    u = zeros(m, 1);
    return;
end

% Construct an LP using scaling parameter to enforce direction preserving.
A = [B -yd];
b = m_higher - B * uMin;
c = [zeros(m, 1); -1];
h = [uMax - uMin; 1];

% To find feasible solution construct problem with appended slack variables.
% A.6.4 Initialization of the Simplex Algorithm of Aircraft control allocation.
sb = 2 * (b > 0) - 1;
Ai = [A diag(sb)];
ci = [zeros(m + 1, 1); ones(n, 1)];
inBi = m + 2:m + n + 1;
ei = true(m + n + 1, 1);
hi = [h; 2 * abs(b)];

% Use Bounded Revised Simplex to find initial basic feasible point of
% original program.
[y1, inB1, e1, itlim, errsimp] = simplxuprevsol_tiebreak(Ai, ci', b, inBi, hi, ei, opts, n, m + n + 1, itlim);

% Check that feasible solution was found.
if itlim <= 0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1 > (m + 1))
    errout = -2;
    disp('No Initial Feasible Solution found'); % m_higher unattainable
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~= 0
    % Construct an incorrect solution to accompany error flags.
    xout = zeros(m + 1, 1);
    indv = inB1 <= (m + 1);
    xout(inB1(indv)) = y1(indv);
    xout(~e1(1:m + 1)) = -xout(~e1(1:m + 1)) + h(~e1(1:m + 1));
else
    % Solve using initial problem from above.
    [y2, inB2, e2, itlim, errsimp] = simplxuprevsol_tiebreak(A, c', b, inB1, h, e1(1:m + 1), opts, n, m + 1, itlim);

    % Construct solution to original LP problem from bounded simplex output.
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

% Transform back to control variables.
u = xout(1:m) + uMin;
lambda = xout(m + 1);

end
