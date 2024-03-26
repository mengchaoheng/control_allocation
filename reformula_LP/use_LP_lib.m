function [u,a] = use_LP_lib(B,v, umin, umax)
% we just convert the origin problem and then use the LP lib.
    % DIR_ALLOC - Direct control allocation.
    %
    %  [u,a] = dir_alloc(B,v,umin,umax)
    %
    % Performs direct control allocation by solving the LP
    %
    %   max a   subj. to  Bu = av
    %   a,u               umin <= u <= umax
    %                        0 <= a
    % If a > 1, set u = u/a.
    %
    % Note: This function has not been optimized for speed.
    %
    %  Inputs:
    %  -------
    % B     control effectiveness matrix (k x m)
    % v     commanded virtual control (k x 1)
    % umin  lower position limits (m x 1)
    % umax  upper position limits (m x 1)
    % 
    %  Outputs:
    %  -------
    % u     optimal control
    % a     scaling factor
    %
    % Reformulate problem to fit linprog format:
    % In matlab, the forms of problem is:
    % min f'x subj. to A*x <= b
    %                  Aeq*b = beq
    %		           lb <= x <= ub
    %% but we use the Standard Forms for Linear Programming Problems
    % min c'x subj. to A*x =b
    %                  0 <= x
    %% so we have to reformula the direction-preserving control allocation
    % problem to:
    % min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
    %                                umin <= u <= umax
    %                                   0 <= a
    % and set x=u-umin, then
    % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
    %                                        x <= umax-umin
    %                                        0 <= x 
    %                                        0 <= a
    % add slack to converted inequalities into equalities of x:
    % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
    %                                     x + s = umax-umin
    %                                         0 <= x 
    %                                         0 <= a
    %                                         0 <= s
    % set X=[x; a; s], that is:
    % min z=[0; -1; 0]'[x; a; s]   s.t.  [B -v 0; I 0 I][x; a; s] = [-B*umin; umax-umin] 
    %                                         0 <= x 
    %                                         0 <= a
    %                                         0 <= s
    % A=[B -v 0; I 0 I]; b=[-B*umin; umax-umin]; c=[0; -1; 0]; X=[x; a; s]

    tol=(1e-8);
    % Number of variables
    [k,m] = size(B);
    % Locate the maximum magnitude element in the desired objective
    [my,~]=max(abs(v));
    
    %Trivial solution, if desired moment is close to zero
    %  May want to adjust the tolerance to improve numerics of later steps
    if (my < tol)    %yd = 0 ==> u=0
        u = (zeros(m,1));
        a=0;
        return;
    end

    %% 1. use the standard form, error!
    A=[B -v zeros(k,m); eye(m) zeros(m,1) eye(m)]; 
    b=[-B*umin; umax-umin]; 
    c=[zeros(m,1); -1; zeros(m,1)];
    Eqin=zeros(k+m,1);
    %% 2. use the general form (LP lib)
    % min c'x subj. to A*x (Eqin) b
    %                     0 <= x
    % [1] N. Ploskas and N. Samaras, Linear Programming Using MATLAB®, vol. 127. in Springer Optimization and Its Applications, vol. 127. Cham: Springer International Publishing, 2017. doi: 10.1007/978-3-319-65919-0.
    % So the param of this lib is: A=[B -v; I 0]; b=[-B*umin; umax-umin]; c=[0; -1]; X=[x; a]
    % A=[B -v; eye(m) zeros(m,1)]; 
    % b=[-B*umin; umax-umin];
    % c=[zeros(m,1); -1];
    % Eqin=[zeros(k,1); -ones(m,1)];

    % [A1, c1, b1, Eqin1, MinMaxLP1] =  general2standard(A, c, b, Eqin, -1);

% [A2, c2, b2, Eqin2, MinMaxLP2] =  standard2canonical(A1, c1, b1, Eqin1, MinMaxLP1);

    %% Use LP lib to analysis.  
    % MinMaxLP=[];
    % c0=[];
    % reinv=[];
    % tole1=[];
    % tole2=[];
    % tole3=[];
    % scalingTechnique=[];
    % pivotingRule=[];
    % basisUpdateMethod=[];
    % 
    % maxIterations=[];
    % tol=[];
    % etaMin=[];
    %% rsa, rdsa and epsa is error allocation
    % Revised Primal Simplex Algorithm
    % [xsol, fval, exitflag, iterations] = rsa(A, c, b, Eqin);
    % % Revised Dual Simplex Algorithm
    % [xsol, fval, exitflag, iterations] = rdsa(A, c, b, Eqin;
    % % Exterior Point Simplex Algorithm
    % [xsol, fval, exitflag, iterations] = epsa(A, c, b, Eqin)

    %% successful !!!
    % % Interior Point Methods
    % [xsol, fval, exitflag, iterations] = ipdipm(A, c, b, Eqin);

     [xsol, fval, exitflag, iterations] = linprogSolver(A, c, b, Eqin,-1,0, 'dual-simplex')
     %%
    % sensitivity Analysis

    % if(exitflag~=1)
    %     u=zeros(m,1);
    %     a=0;
    %     disp('stop!');
    % else
        % x=u-umin, x = xsol(1:m)
        u = xsol(1:m)+umin;
        a = xsol(m+1);
    % end
    % Scale down u if a>1
    if a>1
        u = u/a;
    end
end
