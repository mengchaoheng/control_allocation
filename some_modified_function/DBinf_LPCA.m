function [u, feas, errout,itlim] = DBinf_LPCA(yd,B,wd,up,wu,emax,uMin,uMax,itlim) % note
%   Dual Branch Control Allocation - Linear Program
%      Objective Error Minimization (1-norm)
%      Control Error Minimization (inf-norm)
%
% [u, feas, errout,itlim] = DB_LPCA(yd,B,wd,up,wu,emax,uMin,uMax,itlim);
%
%    Uses a Bounded Revised Simplex solver to solve two linear Programming problems
%  The first branch, feasibility, seeks to minimize the weighted 1-norm of the
%    objective error
%    min J= | diag(wd)*(B*u - yd) |_1  s.t. umin <= u <=umax
%
%   If J = 0 then yd is on the interior of the AMS and a second program seeks to minimize
%   the maximum of a set of weighted absolute errors compared to a "preferred" control solution.
%    min  [ max (i=1:n) wu(i)*abs(u(i)-up(i))]  s.t. umin <= u <= umax and  Bu=yd
%
%   (In the text, the overall structure is presented in A.5.1)
%
%   (The overall structure of this allocator is similar to that presented
%   in Buffington, J. et al. "On-Line System Identification for Aircraft with Distributed
%         Control Effectors", Int. J. Robust Nonlinear Control 9, 1033-1049 (1999), however
%    the control error minimization branch is modified).
%        
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          wd [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution
%          wu [m,1]  = Weighting on control error
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         feas       = Feasible flag
%                        0 = yd is unachievable
%                        1 = Bu=yd
%                        2 = Bu=yd and u = up
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in Feasible branch
%                         >0 = Error in Sufficient branch
%                         -1,1 = Solver error (unbounded solution)
%                         -3,3 = Iteration limit exceeded
%                         -4   = Objective error saturates emax
%         itlim      = Number of iterations remaining after solution found
% 
% Calls:
%        DBcaLP1f_sol  -- Subfunction for feasibility branch
%        DBcaLP1s_sol  -- Subfunction for sufficient  branch
%
% Notes:
%    Feasibility is assessed with a hard-coded tolerance on the cost, J, of 1e-5. This
%        should be reset based on the scaling of the inputs.
%
%    The "preferred" control, up, is used to initialize the feasibility branch. The resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible.Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
%    Error code < 0 implies an error in the first branch and there is no guarantee on
%  the quality of the output solution other than the control limits and 
%     wd'*(B*u-yd) <= wd'*(B*yp-yd)
%   
%    Error code > 0 for errors in second branch and, while it may not be the minimal 
%  control error, the resulting solution should have B*u=yd and can be used.
%  
%
% Modification History
%   2002      Roger Beck  Original (DPcaLP1it.m)
%   8/2014    Roger Beck  Update for use in text
%

%Figure out how big problem is (use standard CA definitions for m & n
[n,m] = size(B);

%Call Feasibility branch
[u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,wd,emax,up,uMin,uMax,n,m,itlim);

feas = 0;
%Check if feasible...if so call sufficiency branch if not exit
if J < 1e-5
    feas = 1;
	[u,Js, errout,itlim] = DBinfcaLP1s_sol(yd,B,wu,u,up,inBout, eout, uMin,uMax,n,m,itlim);   
    if Js < 1e-5
       feas = 2;
    end
end
return;

end % DB_CALPcaLP