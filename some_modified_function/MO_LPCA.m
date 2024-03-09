function [u,errout] = MO_LPCA(yd,B,up,lam, eMax,uMin,uMax,itlim)
% Mixed Optimization (Single Branch) Control Allocation Linear Program
%    Objective Error Minimizing
%    Control Error minimizing
%
% function [u,errout] = MO_LPCA(yd,B,w,up,uMin,uMax,itlim);
%
%    Solves the control allocation problem while seeking to
%  simultaneously minimize the error in the desired objective
%  and the error in the controls using a single linear program.
%
%  Finds the solution that minimizes
%  min |B*u - yd |_1+lambda*|(u-up)|_1
%   such that
%   uMin <= u <= uMax
%
%  The balance between the two objectives is determined by the
%  weight lambda. This weight needs to be chosen small so that
%  the objective error minimization part of the problem dominates
%  the solution in order to ensure that B*u is close to yd for
%  attainable solutions.
%
%  (Section A.5.2 and example A.7 in the text discuss Single Branch
%   optimization routines including this formulation).
%
%   (See Bodson, M., "Evaluation of Optimization Methods for
%          Control Allocation",  AIAA 2001-4223).
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          up[m,1]   = Preferred Control Vector
%          lam       = Control Error Weighting (single parameter)
%          eMax[n,1] = Maximum Objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         -1 = Solver error (unbounded solution)
%                         -3 = Iteration limit exceeded
%                         -4 = Initial feasible solution saturates emax
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   Attainable controls may result in a small error in objective depending on 
%   scaling of control and objective errors and the magnitude of the weight, lambda.
%
%   Solution with errors, errout~=0, only ensures that control limits are respected and
%  that weighted error, |B*u - yd |_1+lambda*|(u-up)|_1, is <= |B*up-yd|_1.
%
%   The "preferred" control, up, is used to initialize the optimization. The resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible. Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
%   If errout ~0 there was a problem in the solution. %
%
% Modification History
%   2002      Roger Beck  Original (MOcaLP1.m(
%   8/2014    Roger Beck  Update for use in text

%Initialize error code to zero
errout = 0;

%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Formulate as an LP problem
A = [eye(n) -eye(n) -B B];
b = B*up-yd;
c = [ones(2*n,1); lam*ones(2*m,1)];
h = [eMax;eMax;uMax-up;up-uMin];

% A feasible initial condition is the up, using the objective error as
%  slack variables.
eyp = B*up-yd;
x0 = [ max(eyp,zeros(n,1));max(-eyp,zeros(n,1));zeros(2*m,1)];

%Find Basis Variables for initial solution
%   If preferred control has zero objective error in an axis, identify
%     additional basic variables that are zero.
%   This is unlikely with floating point data and limited precision
%     --could also handle by biasing zero terms in eyp by eps.
%  
%
indn = 1:n;
numzer = length(find(eyp==0));
inBi = [indn(eyp>0) n+indn(eyp<0) (2*n):( (2*n)-1+numzer )];
e = true(2*(n+m),1);

%Solve using Bounded Revised Simplex
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inBi,h,e,n,2*(m+n),itlim);

   %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(2*(m+n),1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
       
    %Check if solution contains error terms that are limited at their upper limit
tmp = ~e2;
tmp(inB2) = false;
if any(tmp(1:2*n))
   errout = -4;
   disp('Output objective error at bounds');
end

    if itlim<=0
        errout = -3;
        disp('Too Many Iterations Finding Final Solution');
    end
    if errsimp
        errout = -1;
        disp('Solver error');
    end

%%Transform Solution Back Into control variables
% Note that x(2*n+1:2*n+m) are the + differences from the preferred control
% solution and x(2*n+m+1:2*(n+m)) are the - differences.

u = xout(2*n+1:2*n+m)-xout(2*n+m+1:2*(n+m))+up;
return;
end
