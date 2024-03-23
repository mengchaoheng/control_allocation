function [u,a] = dir_alloc_linprog_re_bound(B,v, umin, umax, lam)
  
    % DIR_ALLOC - Direct control allocation.
    %
    %  [u,a] = dir_alloc(B,v,umin,umax)
    %
    % Performs direct control allocation by solving the LP
    %
    %   max a   subj. to  Bu = av
    %   a,u               umin <= u <= umax
    %                       0 <= a <= lam
    % lam >= 1,if lam = 1, it will at limit, but the allocation error is 0.
    % when lam > 1, if a > 1, set u = u/a. 
    % lam can not be set to Inf, matlab will stop. It can be fix by
    % reformula the problem, remove the bound of a, become dir_alloc_linprog_re.
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
    % min f'x subj. to A*x <=b
    %                  Aeq*b = beq
    %		           lb <= x <= ub
    %% but we use the Standard Forms for Linear Programming Problems
    % min c'x subj. to A*x =b
    %                  0 <= x
    %% so we have to reformula the direction-preserving control allocation
    % problem to:
    % min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
    %                                umin <= u <= umax
    %                                   0 <= a <= lam
    % and set x=u-umin, then
    % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
    %                                0 <= x <= umax-umin
    %                                0 <= a <= lam
    % add slack to converted inequalities into equalities of x and a:
    % min z=[0; -1]'[x; a]   s.t.     [B -v][x; a] = -B*umin
    %                            [x; a] + [sx; sa] = [umax-umin; lam]
    %                                            0 <= x 
    %                                            0 <= a
    %                                            0 <= sx
    %                                            0 <= sa
    % set X=[x; a; sx; sa], that is:
    % min z=[0; -1; 0; 0]'[x; a; sx; sa]   s.t.  [B -v 0 0; I 0 I 0][x; a; sx; sa] = [-B*umin; umax-umin; lam] 
    %                                            0 <= x 
    %                                            0 <= a
    %                                            0 <= sx
    %                                            0 <= sa
    % A=[B -v 0 0; I 0 I 0]; b=[-B*umin; umax-umin; lam]; c=[0; -1; 0; 0];
    
    % Number of variables
    [k,m] = size(B);
    A=[B -v zeros(k,m) zeros(k,1); eye(m) zeros(m,1) eye(m) zeros(m,1); zeros(1,m) 1 zeros(1,m) 1]; 
    b=[-B*umin; umax-umin; lam]; % a的上界lam=1，但此处设为可调。
    c=[zeros(m,1); -1; zeros(m,1); 0];
    
    % Solve linear program
    % Algorithm = 'dual-simplex'（默认值）
    % 
    %              'interior-point-legacy'
    % 
    %              'interior-point'
    % Display = 'final'（默认值）仅显示最终输出。
    % 
    %             'off' 或 'none' 不显示输出。
    % 
    %             'iter' 在每次迭代时显示输出。

    options = optimset('Display', 'off');
    [X,fval,exitflag,output,lambda]= linprog(c,[],[],A,b,zeros(2*m+2,1),[],options);
    if(exitflag~=1)
        u=zeros(m,1);
        a=0;
        disp('stop!');
    else
        % x=u-umin, x = X(1:m)
        u = X(1:m)+umin;
        a = X(m+1);
    end
    % Scale down u if a>1
    if a>1
        u = u/a;
    end
end
