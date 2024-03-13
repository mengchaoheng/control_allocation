%% Iterate through simplex algorithm main loop
function [x,z,iters]=Simplex_loop_mch(B, A, b, c, z)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=、 >=、 <=) b
%   x 
    %% Initialization
    tol=1e8;    
    iters=0;
    [m,n] = size(A);
    while ~all(c>=-1/tol) % 0 is 1/tol               % 3.~isempty(c(c(N)<0))
    e = find(c < -1/tol, 1, 'first'); % 进基变量索引    % 4. e = N(find(c(N)<0,1))
    a_ie=A(:,e);
    ip=a_ie>(1/tol);
    delta=tol*ones(m,1);
    if ~isempty(ip)
        delta(ip)=b(ip)./a_ie(ip);
    end
    [~,L]=min(delta);%选择离基 (离基在B数组中的行索引)
%         li = B(L);    % 离基变量索引               
    if delta(L) >= tol                         
        disp('System is unbounded!');
        return;  % 11.
    else
        [B,A,b,c,z] = pivot_mch(B,A,b,c,z,L,e); %换基，即进行初等行变换
    end
    iters = iters + 1;
       
    end
    x = zeros(n,1);                               % [13,16].
    x(B) = b; 
end
