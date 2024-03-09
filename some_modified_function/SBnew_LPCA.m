function [u,errout] = SBnew_LPCA(yd,B,w,up,uMin,uMax) % note
% Single Branch Control Allocation Linear Program
%    Direction Preserving
%    Control Error minimizing
%
% function [u,errout] = SB_LPCA(yd,B,w,up,uMin,uMax,itlim);
%
%    Solves the control allocation problem while seeking to
%  simultaneously preserve the direction for unattainable
%  objectives and minimizing the control error for attainable
%  commands using a single linear program.
%
%  Finds the solution that minimizes
%  min -lambda + |diag(w)*(u-up)|_1
%   such that
%  B*u = lambda*yd
%   uMin <= u <= uMax
%      0 <= lambda <= 1
%
%  The balance between the two objectives is determined by the
%  weight vector, w. For unattainable moments the constraints
%  ensure that the direction of the command is maintained, 
%  however, the weights should be small compared to the
%  relative scaling of the units on the control vector and the
%  objective vector, so that the error term doesn't pull the
%  optimum solution away from B*u =yd.
%
%  (Section A.5.2 and example A.6 in the text discuss Single Branch
%  optimization routines of including this formulation).
% 
%   (See Buffington, J. "Tailess Aircraft Control Allocation",
%      AIAA-97-3695 for a similar approach that seeks to minimize
%      control usage (i.e. up = 0) and partitions lambda to prioritize
%      command components)
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [m,1]   = Control Error Weighting
%          up[m,1]   = Preferred Control Vector
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                       1  linprog converged to a solution X.
%                       0  Maximum number of iterations reached.
%                      -2  No feasible point found.
%                      -3  Problem is unbounded.
%                      -4  NaN value encountered during execution of algorithm.
%                      -5  Both primal and dual problems are infeasible.
%                      -7  Magnitude of search direction became too small; no further
%                           progress can be made. The problem is ill-posed or badly
%                           conditioned.
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         linprog
%
% Notes:
%   Attainable controls may result in a small error in objective depending on 
%   scaling of control and objective errors and the magnitude of the weights.
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution, result has B*u in the right direction
% and magnitude <= yd.
%
% Modification History
%   2002      Roger Beck  Original (SBcaLP2)
%   8/2014    Roger Beck  Update for use in text


%Initialize error code to zero
errout = 0;

%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

A=[B           -B         -yd      zeros(n,m) zeros(n,m) zeros(n,1);
   eye(m)    zeros(m,m)  zeros(m,1) eye(m)    zeros(m,m) zeros(m,1);
   zeros(m,m) eye(m)     zeros(m,1) zeros(m,m) eye(m)    zeros(m,1);
   zeros(1,m) zeros(1,m) 1          zeros(1,m) zeros(1,m) 1];
b=[[0;0;0]-B*up;uMax-up;up-uMin;1];
c=[w;w;-1;zeros(m,1); zeros(m,1); 0];
%---------------------ok----------------------
[x,~,exitflag,~,~] = linprog(c',[],[],A,b,zeros(size(A,2),1),ones(size(A,2),1)*1000);
%--------------------
 u=x(1:m)-x(m+1:2*m);
 errout=exitflag
end