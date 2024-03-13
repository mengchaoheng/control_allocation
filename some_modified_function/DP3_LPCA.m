function [u, errout,rho] = DP3_LPCA(yd,B,uMin,uMax,itlim)
% Direction Preserving Control Allocation Linear Program
%
% function [u, errout] = DP_LPCA(yd,B,uMin,uMax,itlim);
%
%    Solves the control allocation problem while preserving the
%  objective direction for unattainable commands. The solution
%  is found by solving the problem,
%    min -lambda,
%    s.t. B*u = lambda*yd, uMin<=u<=uMax, 0<= lambda <=1
%
%  For yd outside the AMS, the solution returned is that the
%  maximum in the direction of yd.
%
%  For yd strictly inside the AMS, the solution achieves
%  Bu=yd and m-n controls will be at their limits; but there
%  is no explicit preference to which solution will be 
%  returned. This limits the usefulness of this routine as
%  a practical allocator unless preferences for attainable solutions
%  are handled externally.
%
%  (For derivation of a similar formulation see A.1.2 and A.2.3 in the
%  text)
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   If errout ~0 there was a problem in the solution. %
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution--B*u is in the correct direction and has
%  magnitude < yd, but B*u may not equal yd (for yd attainable)
%   or be maximized (for yd unattainable)
%
% Modification History
%   2002      Roger Beck  Original (DPcaLP8.m)
%   8/2014    Roger Beck  Update for use in text


%Initialize error code to zero
tol = 1e-1;
%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);
%Check to see if yd == 0
%  May want to adjust the tolerance to improve numerics of later steps
% if (all(abs(yd) < tol))    %yd = 0 ==> u=0
if norm(yd) < tol
    errout = -1;
    u = pinv(B)*yd;% zeros(m,1);
    rho=2; % 大于1
    return;
end
%Construct an LP using scaling parameter to enforce direction preserving
% A = [B -yd];
% b = -B*uMin;
% c = [zeros(m,1);-1];
% h = [uMax-uMin; 1e4];
% %---------------------ok----------------------
% [x,~,exitflag,~,~] = linprog(c',[],[],A,b,zeros(size(A,2),1),h);
% %--------------------
% u = x(1:m)+uMin;
% errout=exitflag;
% rho = x(m+1);


A=[B -yd zeros(n,m) zeros(n,1);
   eye(m) zeros(m,1) eye(m) zeros(m,1);
   zeros(1,m) 1 zeros(1,m) 1];
b=[-B*uMin;uMax-uMin; 1];
c = [zeros(m,1);-1;zeros(m,1);0];
% 换第三方解法

[x,~,exitflag,~,~] = linprog(c',[],[],A,b,zeros(size(A,2),1),[]);

u = x(1:m)+uMin;
errout=exitflag;
rho = x(m+1);
% if rho > 1
%     u = u/rho;
% end
end