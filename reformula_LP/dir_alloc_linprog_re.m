function [u,a] = dir_alloc_linprog_re(B,v, umin, umax)
  
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
    %                                   0 <= a
    % and set x=u-umin, then
    % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
    %                                0 <= x <= umax-umin
    %                                0 <= a
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
    % A=[B -v 0; I 0 I]; b=[-B*umin; umax-umin]; c=[0; -1; 0];
    
    
    % Number of variables
    [k,m] = size(B);
    A=[B -v zeros(k,m); eye(m) zeros(m,1) eye(m)]; 
    b=[-B*umin; umax-umin]; 
    c=[zeros(m,1); -1; zeros(m,1)];
    
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
    [X,fval,exitflag,output,lambda]= linprog(c,[],[],A,b,zeros(2*m+1,1),[],options);
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
