function [u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,w,emax,up,uMin,uMax,n,m,itlim)
%   Dual Branch Control Allocation--Feasibility Branch
%      Objective Error Minimization
%
% [u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,w,emax,up,uMin,uMax,n,m,itlim);
%
%    Determines if a feasible control solution can be found minimizing the error between
%  the desired and obtained objective (i.e. is yd in the AMS?). If yd is unattainable
%  the returned solution minimizes the weighted 1-norm of the objective error.
%  The Bounded Revised Simplex solver is called to minimize
%    min J= | diag(wd)*(B*u - yd) |_1  s.t. umin <= u <=umax
%   (See A.1.2 and Example A.2 in the text for a discussion of a similar formulation).
%
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution (used for initialization)
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          n         = Number of objectives
%          m         = Number of controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         J          = Cost=weighted 1-norm of objective error wd'*abs(Bu-yd)
%         eout       = Bound flag for controls (true at lower bound, false at upper bound)
%         errout     = Error Status code
%                         0 = found solution
%                        -3 = Iteration limit exceeded in branch
%                        -1 = Solver error
%                        -4 = Objective error for solution has at least one component
%                               at emax.
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    The "preferred" control, up, is used to initialize the solver and the resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible.Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update
%

%Initialize error code to zero
errout = 0;

%Formulate as an LP problem
A = [eye(n) -eye(n) -B];
b = B*uMin-yd;
c = [w;w;zeros(m,1)];
h = [emax; emax; uMax-uMin];


% A feasible initial condition is the up, using the objective error as
%  slack variables.
eyd = B*up-yd;
x0 = [ max(eyd,zeros(n,1));max(-eyd,zeros(n,1));zeros(m,1)];

%Find Basis Variables for initial solution
%   If preferred control has zero objective error in an axis, identify
%     additional basic variables that are zero.
%   This is unlikely with floating point data and limited precision
%     --could also handle by biasing zero terms in eyp by eps.
%
%
indn = 1:n;
numzer = length(find(eyd==0));
inBi = [indn(eyd>0) n+indn(eyd<0) (2*n):( (2*n)-1+numzer )];

e = true(2*n+m,1);

%Solve using Bounded Revised Simplex
[y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inBi,h,e,n,2*n+m,itlim);

%Construct solution to original LP problem from bounded simplex output
%  Set non-basic variables to 0 or h based on e2
%  Set basic variables to y2 or h-y2.
xout = zeros(2*n+m,1);
xout(inB2) = y2;
xout(~e2) = -xout(~e2)+h(~e2);


if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding Final Solution');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

%Check if solution contains error terms that are limited at their upper limit
tmp = ~e2;
tmp(inB2) = false;
if any(tmp(1:2*n))
   err = -4;
   disp('Output objective error at bounds');
end
%Compute cost (objective error)-- used to determine feasibility in
%subsequent program
J = c'*xout;

%Convert solution back to control variable
u = xout(2*n+1:2*n+m)+uMin;

%Output controls in basis for subsequent program
%  If solution is feasible, then any error terms in basis are 0
%   and can be replaced with one of the controls at zero

inBout = inB2(inB2>(2*n))-2*n;
if length(inBout) < n  %If there are too few controls in basis
    cind = 1:m;
    cvec = setdiff(cind,inBout);
    inBout = [inBout cvec(1:(n-length(inBout)))];
end

eout = e2(2*n+1:2*n+m);

end %DBcaLP1f_sol