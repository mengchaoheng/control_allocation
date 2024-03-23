function [u,a] = dir_alloc_simplex_re(B,v, umin, umax)
  
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
    % Reformulate problem to fit LP lib format:
    % In LP lib, the forms of problem is:
    % max f'x      s.t.  A*x <=b
    %                    Aeq*x = beq
    %		             0 <= x 
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
    % Aeq=[B -v]; beq=[-B*umin]; A=[I 0]; b=[umax-umin]; c=[0; -1];  X=[x; a]
    
    % Number of variables
    [~,m] = size(B);
    Aeq=[B -v]; 
    beq=-B*umin;
    A=[eye(m) zeros(m,1)]; 
    b=umax-umin; 
    c=[zeros(m,1); -1];
    %% use linprog
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
    [X,fval,exitflag,output,lambda]= linprog(c,A,b,Aeq,beq,zeros(m+1,1),[],options);
    if(exitflag~=1)
        u=zeros(m,1);
        a=0;
        disp('stop!');
    else
        % x=u-umin, x = X(1:m)
        u = X(1:m)+umin;
        a = X(m+1);
    end
    %% ToDo: use LP lib (git@github.com-mch:mengchaoheng/control_allocation.git)
    
    u = X(1:m)+umin;
    a = X(m+1);
    %% Scale down u if a>1
    if a>1
        u = u/a;
    end
end
