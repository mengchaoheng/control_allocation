function [u] = LPwraparm(IN_MAT)
% Make single input, single output version of Linear Programming for use in
% Simulink via the MATLAB Fcn block
% IN_MAT = [B     d         
%           umin' 0        
%           umax' 0]
%
% 20140905  Created version to use Roger Beck's DB_LPCA program
% 20151206  Updated with Roger Beck's latest version of code and added 
%           LPmethod to select the various Linear Programming algorithms 
% Get sizes
% If matrices too small, set contols to zero and return
if  norm(IN_MAT)<1e-16
    u=zeros(4,1);
    return
end
% Partition input matrix into component matrices
B=IN_MAT(1:3,1:4);
yd=IN_MAT(1:3,5);
uMin=IN_MAT(4,1:4)';
uMax=IN_MAT(5,1:4)';
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
%          lam       = Control Error Weighting (single parameter)
%          eMax[n,1] = Maximum Objective error
%          w [m,1]   = Control Error Weighting

% change variable names to be consistent with Roger's documentation
% Some variables not specifically defined for generic case
% Values below are suggestions for now
% Users may add code here to specify alternative values
itlim=50;
[u_act, ~] = DParm_LPCA(yd,B,uMin,uMax,itlim,3,4);%DP_LPCA(yd+ye-Bu0,B,uMin,uMax,itlim)
        % Direction Preserving Control Allocation Linear Program    
u=u_act;
end
function [y0, inB, e,itlim,errout] = simpl(A,ct,b,inBx,h,e,varargin)
%  Bounded Revised Simplex
%
%function [yout, inBout,bout, itout,errout] = simplxuprevsol(A,ct,b,inB,inD,h,e,m,n,itlim)
%
%   Solves the linear program:
%          minimize c'y 
%          subject to 
%          Ay = b
%          0<= y <= h
%
%  Inputs: 
%          A [m,n]   = lhs Matrix of equaltity constraints
%          ct [1,n]  = transpose of cost vector
%          b [m,1]   = rhs vector for equality constraint
%-----------------------------------------------------------------------
%          inB [m]   = Vector of indices of unknowns in the initial basic set
%          inD [n-m] = Vector of indices of unknowns not in the initial basic set
% from inBx [n], inBx is a logical vector that non-zeros element indices
% are inB, so as inD
%--------------------------------------------------------------------------------
%          h[n,1]    = Upper Bound for unknowns
%          e[n,1]    = Sign for unknown variables (+ lower bound, - upper bound)
%  Optional inputs:
%          m,n       = number of constraints, unknowns (Opposite standard
%                      CA convention
%          itlim     = Upper bound on the allowed iterations\
%
% Outputs:
%         yout[n,1]  = Optimal output variable
%         inBout     = indices of Basic vectors in output
%         eout       = sign associate with output unknowns
%         itout      = number of iterations remaining out of itlim
%         errout     = Flag (=true) if unbounded is set
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update for use
%   9/2014    Roger Beck  Added anti-cycling rule

%Optional Inputs
switch  length(varargin)
    case 0
    itlim = inf;
    [m,n] = size(A);
    case 1
    itlim = varargin{1};
    [m,n] = size(A);
    case 2
    itlim = inf;
    m = varargin{1};
    n = varargin{2};
    case 3
    itlim =varargin{3};
    m = varargin{1};
    n = varargin{2};
 end    	
    	
%Tolerance for unknown == 0
tol = 1e-10;

%Index list for non-basic variables
nind = 1:(n-m);

%Partition A
inB=find(inBx);
inD=find(~inBx);
%-----------------------
% inD = setdiff(1:n, inB);
%-----------------------

%Adjust signs problem if variables are initialized at upper
% bounds.
A(:,~e) = -A(:,~e);
ct(~e) = -ct(~e);
b = b + A(:,~e)*h(~e);

y0 = A(:,inB)\b;  %Initial Solution

%Initialize Loop Termination Conditions
done = false;
unbounded = false;

%Main Simplex loop
while (~done  || ~unbounded ) && (itlim > 0)
    itlim = itlim-1;

    %Calculate transpose of relative cost vector based on current basis
    lamt = ct(inB)/A(:,inB);
    rdt = ct(inD)-lamt*A(:,inD);
    %Find minimum relative cost
    [minr, qind] = min(rdt);
    if minr >=0  % If all relative costs are positive then the solution is optimal
        done = true;
        break;
    end
    qel = inD(qind);  % Unknown to Enter the basis minimizes relative cost
    yq = A(:,inB)\A(:,qel); %Vector to enter in terms of the current Basis vector
    
    if all(abs(yq)<=tol)
      unbounded = true;
      disp(' Solution is unbounded');  % Check this condition
      break
    end

    %Compute ratio how much each current basic variable will have to move for the entering
    % variable.

    rat = y0./yq; 
    
    % If yq < 0 then increasing variable when it leaves the basis will minimize cost
    hinB = h(inB);
    indm = yq<0;
    rat(indm) = rat(indm) - hinB(indm)./yq(indm);
    % If an element yq ~=0 then it doesn't change for the entering variable and shouldn't
    %  be chosen
    indz = abs(yq)<=tol;
    rat(indz) = inf;

    % Variable to exit is moving to its minimum value
    [minrat, p] = min(rat);

   % If the minimum ratio is zero, then the solution is degenerate and the entering
   %   variable will not change the basis---invoke Bland's selection rule to avoid
   %   cycling.
    if (abs(minrat) <= tol)
       % Find negative relative cost
       indm = nind(rdt<0); %Note that since minr <0 indm is not empty   
       qind = indm(1);
       qel = inD(qind);  % Unknown to Enter the basis is first indexed to avoid cycling
       yq = A(:,inB)\A(:,qel); %Vector to enter in terms of the current Basis vector
       if all(abs(yq)<=tol)
           unbounded = true;
           disp(' Solution is unbounded');  % Check this condition
           break
       end
       % Recompute rations and determine variable to leave
       rat = y0./yq; 
        % If yq < 0 then increasing variable when it leaves the basis will minimize cost
        hinB = h(inB);
        indm = yq<0;
        rat(indm) = rat(indm) - hinB(indm)./yq(indm);
        % If an element yq ~=0 then it doesn't change for the entering variable and shouldn't
        %  be chosen
        indz = abs(yq)<=tol;
        rat(indz) = inf;

        % Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        [minrat, p] = min(rat);
    end

  % Maintain the bounded simplex as only having lower bounds by recasting 
  % any variable that needs to move to its opposite bound.
    if (minrat >= h(qel))
           %Case 1: Entering variable goes to opposite bound and current basis is maintained
            e(qel) = ~e(qel);
            A(:,qel) = -A(:,qel);
             b = b + A(:,qel)*h(qel);
             ct(qel) = -ct(qel);
    elseif yq(p) > 0
           %Case 2: Leaving variable returns to lower bound (0)	
           pel = inB(p);
           inB(p)= qel;
           inD(qind)= pel;
     else
           %Case 2: Leaving variable moves to upper bound	
            pel = inB(p);
            e(pel)=~e(pel);
            A(:,pel) = -A(:,pel);
            inB(p)= qel;
            inD(qind)= pel;
            ct(pel) = -ct(pel);
            b = b + A(:,pel)*h(pel);
     end
        
    y0 = A(:,inB)\b; % Compute new Basic solution;
end
errout = unbounded;     
end
function [u, errout] = DParm_LPCA(yd,B,uMin,uMax,itlim, n,m)
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
errout = int32(0);

tol = 1e-10;
%Figure out how big the problem is (use standard CA definitions for m & n)
% [n,m] = size(B);

%Check to see if yd == 0
%  May want to adjust the tolerance to improve numerics of later steps
if (all(abs(yd) < tol))    %yd = 0 ==> u=0
    errout = int32(-1);
    u = zeros(m,1);
    return;
end

%Construct an LP using scaling parameter to enforce direction preserving
A = [B -yd];
b = -B*uMin;
c = [zeros(m,1);-1];
h = [uMax-uMin; 1];


%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(m+1,1);ones(n,1)];
% inBi = m+2:m+n+1;
inBi=false(1,m+n+1);
inBi(m+2:m+n+1) =1;
ei = true(m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point of
%original program
% [y1, inB1, e1,itlim, errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,m+n+1,itlim);
[y1, inB1, e1,itlim, errsimp] = simpl(Ai,ci',b,inBi,hi,ei,n,m+n+1,itlim);
%Check that Feasible Solution was found
if itlim<=int32(0)
    errout = int32(-3);
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(m+1))
    errout =int32(-2);
    disp('No Initial Feasible Solution found');
end
	if errsimp
	    errout = int32(-1);
		disp('Solver error');
	end

if errout ~=int32(0)  % Construct an incorrect solution to accompany error flags
        xout = zeros(m+1,1);
        indv = inB1<=(m+1);
        xout(inB1(indv)) = y1(indv);
        xout(~e1(1:m+1)) = -xout(~e1(1:m+1))+h(~e1(1:m+1));
else  % No Error continue to solve problem
    
    
    %Solve using initial problem from above
    inBx=false(1,m+1);
    inBx(inB1) =1;
    [y2, inB2, e2,itlim,errsimp] = simpl(A ,c',b,inBx,h,e1(1:m+1),n,m+1,itlim);
%     [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:m+1),n,m+1,itlim);
    
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(m+1,1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
    
    if itlim<=int32(0)
        errout = int32(3);
        disp('Too Many Iterations Finding Final Solution');
    end
	if errsimp
	    errout = int32(1);
		disp('Solver error');
	end
    
end


%Transform back to control variables
u = xout(1:m)+uMin;
rho = xout(m+1);
if rho > 1
    u = u/rho;
end
return;
end