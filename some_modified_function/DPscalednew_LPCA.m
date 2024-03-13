function [u,itlim,errout,rho] = DPscalednew_LPCA(yd,B,uMin,uMax,itlim)
% Direction Preserving Control Allocation Linear Program
%     Reduced formulation (Solution Scaled from Boundary)
%
% function [u,itlim,errout] = DPscaled_LPCA(yd,B,uMin,uMax,itlim);
%
%    Solves the control allocation problem while preserving the
%  objective direction for unattainable commands. The reduced
%  dimension of the linear program passed to the Bounded Revised
%  Simplex solver is formed by forcing the solution to be on the
%  boundary of the AMS and eliminating the highest magnitude
%  objective by solving the other constraints in terms of it.
%
%  For yd outside the AMS, the solution returned is that the
%  maximum in the direction of yd
%    B*u= lamda*yd
%    max lamda s.t. uMin <= u <= uMax
%
%  Reducing the degrees of freedom elminates the problems of redundant
%  solutions for attainable objectives. If the desired objective is on the
%  interior of the AMS the solution is scaled from the solution on the
%  boundary, yielding the same controls as the Direct Allocation solution.
%  
%  (In the text this solution is discussed in section A.5.3)
%
%   (See Bodson, M., "Evaluation of Optimization Methods for
%          Control Allocation",  AIAA 2001-4223).
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
%    If yd is close to zero, u = 0;
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution.
%
% Modification History
%   2002      Roger Beck  Original ( DPcaLP2.m)
%   8/2014    Roger Beck  Update


%Initialize error code to zero
errout = 0;

tol = 1e-10;
%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

% Locate the maximum magnitude element in the desired objective
[my,iy]=max(abs(yd));

%Trivial solution, if desired moment is close to zero
%  May want to adjust the tolerance to improve numerics of later steps
if (my < tol)    %yd = 0 ==> u=0
    errout = -1;  %Set flag to let caller know that it wasn't solved
    u = zeros(m,1);
    rho=1000;
    return;
end

%Transform Problem by Reordering Objectives with maximum first
Bt = B([iy setdiff(1:n,iy)],:);
ydt = yd([iy setdiff(1:n,iy)]);
ydt(2:3) = ydt([3 2]);Bt([2 3],:) = Bt([3 2],:);
%%Convert into a LP problem
M = [ydt(2:n) -ydt(1)*eye(n-1)];
A = M*Bt;
b = -A*uMin;
c = -Bt'*ydt;
h = uMax-uMin;


% 解方程
[xout,~,exitflag,~,~] = linprog(c',[],[],A,b,zeros(size(A,2),1),h);



%Transform Solution Back Into control variables
% u(i) = x(i)+umin(i) if e(i)
u = xout+uMin;
%Rescale controls so solution is not on boundary of Omega.
rho = ydt'*Bt*u/(ydt'*ydt);
if rho > 1
    u = u/rho;
end


return;
end