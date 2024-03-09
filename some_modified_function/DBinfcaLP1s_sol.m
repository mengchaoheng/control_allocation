function [u,J,errout,itlim] = DBinfcaLP1s_sol(yd,B,w,u0,up,inBi, ei, uMin,uMax,n,m,itlim)
%   Dual Branch Control Allocation--Sufficient Branch
%      Control Error Minimization
%
% function [u,J,errout,itlim] = DBcaLP1s_sol(yd,B,w,u0,up,inBi, ei, uMin,uMax,n,m,itlim);
%
%   Assumes that the desired objective, yd, is attainable and seeks to minimize the maximum
%   individual component of the absolute error between the controls and a "preferred" control solution.
%   (See A.4.2 in the text for a discussion of a similar formulation).
%  The Bounded Revised Simplex solver is called to minimize
%    min J= | diag(wu)*(u - up) |_inf  s.t. umin <= u <=umax
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [n,1]  = Weighting on Control error
%          u0 [m,1]  = Initial basic solution that attains Bu=yd
%          up [m,1]  = Preferred control solution (used for initialization)
%          inBi[n,1] = List of controls in the basis set(i.e. not at limits)
%          ei[n,1]   = Controls not at upper bounds (ei(i) = false if ui(i) = uMax(i))
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          n         = Number of objectives
%          m         = Number of controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         J          = Cost=weighted inf-norm of control error max w(i)*abs(u(i)-up(i))
%         errout     = Error Status code
%                         0 = found solution
%                         3 = Iteration limit exceeded in  branch
%                         1 = Solver error
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    The input solution u0 is assumed to be a basic feasible solution (i.e. yd is attainable
%    and m-n controls are at a limit).
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update 
%   8/2015    Roger Beck  Update from 1-norm to inf-norm 
%


%Figure out how big problem is (use standard CA definitions for m & n
%[n,m] = size(B);

%Initialize Error Code
errout = 0;

%Formulate as an LP problem
A = [B -B zeros(n,2*m+1); ...
     diag(w) zeros(m) eye(m) zeros(m) -ones(m,1) ; ...
     zeros(m) diag(w) zeros(m) eye(m) -ones(m,1); ...
     ];
b = [yd-B*up; zeros(2*m,1)];
c = [zeros(4*m,1); 1];
h = [uMax-up;up-uMin; uMax-up;up-uMin; max(max(uMax-up),max(up-uMin))];

%Any feasible solution that satisfies Bu = yd will satisfy constraint
eup = u0-up;
[us,is] = max([eup; -eup]);
x01 = [max(eup,zeros(m,1));max(-eup,zeros(m,1))];
x02 = us-x01;
x0 = [x01;x02;us];

%Use identified basis from feasible problem--note that limited variables
%above are not in the basis.

inBt = 2*m + (1:2*m);
inBt(is) = [];
inB = [inBi(eup(inBi)>=0) inBi(eup(inBi)<0)+m inBt 4*m+1];


%Determine which non-basic variables are at upper bound
%
e = (x0 < h/2);
e(inB) = true;

%Solve using Bounded Revised Simplex
%e = true(2*m,1);
%e(evec) = false;

[y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB,h,e,n,2*m,itlim);


%Construct solution to original LP problem from bounded simplex output
%  Set non-basic variables to 0 or h based on e2
%  Set basic variables to y2 or h-y2.
xout = zeros(4*m+1,1);
xout(inB2) = y2;
xout(~e2) = -xout(~e2)+h(~e2);

if itlim<=0
    errout = 3;
    disp('Too Many Iterations Finding Final Solution');
end
if errsimp
    errout = 1;
    disp('Solver error');
end
%Compute cost output
J = c'*xout;
%Convert solution back to control variable
u = xout(1:m)-xout(m+1:2*m)+up;

end %DBinfcaLP1s_sol