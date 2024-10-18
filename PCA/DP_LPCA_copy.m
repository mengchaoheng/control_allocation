
function [u, errout,lambda] = DP_LPCA_copy(m_higher,yd,B,uMin,uMax,itlim)
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
errout = 0;
lambda=0;
%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Check to see if yd == 0
%  May want to adjust the tolerance to improve numerics of later steps
if (all(abs(yd) < eps))    %yd = 0 ==> u=0
    errout = -1;
    u = zeros(m,1);
    return;
end

%Construct an LP using scaling parameter to enforce direction preserving
A = [B -yd];
b = m_higher-B*uMin;
c = [zeros(m,1);-1];
h = [uMax-uMin; 1];



%To find Feasible solution construct problem with appended slack variables
% A.6.4 Initialization of the Simplex Algorithm of <Aircraft control allocation>
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(m+1,1);ones(n,1)];
inBi = m+2:m+n+1;
ei = true(m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point of
%original program
[y1, inB1, e1,itlim, errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,m+n+1,itlim);

%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(m+1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
	if errsimp
	    errout = -1;
		disp('Solver error');
    end

if errout ~=0  % Construct an incorrect solution to accompany error flags
    xout = zeros(m+1,1);
    indv = inB1<=(m+1);
    xout(inB1(indv)) = y1(indv);
    xout(~e1(1:m+1)) = -xout(~e1(1:m+1))+h(~e1(1:m+1));
    
else  % No Error continue to solve problem
    
    % or set inB1=[1 2 3] by observe
    %Solve using initial problem from above
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:m+1),n,m+1,itlim);
    
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(m+1,1);
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


%Transform back to control variables
u = xout(1:m)+uMin;
lambda= xout(m+1);
return;
end


