function [u,errout] = SBprio_LPCA(yd,ye,B,w,up,uMin,uMax,itlim) % note
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

%Tolerance for unknown == 0
tol = 1e-10;

%Initialize error code to zero
errout = 0;

%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Formulate as an LP problem
A = [B -B -yd];
b = -B*up+ye;
c = [w; w; -1];
h = [uMax-up; up-uMin;1];

%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(2*m+1,1);ones(n,1)];
% inBi = [2*m+2:2*m+n+1];
inBi=false(1,2*m+n+1);
inBi(2*m+2:2*m+n+1) =1;
ei = true(2*m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point
% [y1, inB1, e1,itlim,errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,2*m+n+1,itlim);
[y1, inB1, e1,itlim,errsimp] = simpl(Ai,ci',b,inBi,hi,ei,n,2*m+n+1,itlim);
%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(2*m+1))
% if any(find(inB1)>(2*m+1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~=0  % Construct an incorrect solution to accompany error flags
    if all( abs(ye)<=tol )
        xout = zeros(2*m+1,1);
        indv = inB1<=(2*m+1);
        xout(inB1(indv)) = y1(indv);
        xout(~e1(1:2*m+1)) = -xout(~e1(1:2*m+1))+h(~e1(1:2*m+1));
    else
%         [u,~,~] = DPscaled_LPCA(ye,B,uMin,uMax,itlim);
        [u,~] = SBprio_LPCA(ye,[0;0;0],B,w,up,uMin,uMax,50);
        return;
    end
else  % No Error continue to solve problem
    
        
    %Solve using initial problem from above
%     [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:(2*m+1)),n,(2*m+1),itlim);
    inBi=false(1,2*m+1);
    inBi(inB1) =1;
    [y2, inB2, e2,itlim,errsimp] = simpl(A,c',b,inBi,h,e1(1:(2*m+1)),n,2*m+1,itlim);
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(2*m+1,1);
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

    
end
    
%Transform Solution Back Into control variables
% Note that x(1:m) are the + differences from the preferred control
% solution and x(m+1:m) are the - differences.

u = xout(1:m)-xout(m+1:2*m)+up;
return;
end