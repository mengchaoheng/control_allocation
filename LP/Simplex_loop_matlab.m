%% Iterate through simplex algorithm main loop
function [x,z]=Simplex_loop_matlab(B, A, b, C, z)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=C*x   subj. to  A*x (=、 >=、 <=) b
%   x 
    %% Initialization
    
    iters=0;
    [m,n] = size(A);
    while ~all(C>=0)                      % 3.~isempty(C(C(N)<0))
        e = find(C < 0, 1, 'first'); % 进基变量索引    % 4. e = N(find(C(N)<0,1))
        a_ie=A(:,e);
        ip=a_ie>0;
        delta=Inf*ones(m,1);
        if ~isempty(ip)
            delta(ip)=b(ip)./a_ie(ip);
        end
        [~,L]=min(delta);%选择离基 行索引
%         li = B(L);    % 离基变量索引               
        if delta(L) == Inf                         
            disp('System is unbounded!');
            return;  % 11.
        else
            [B,A,b,C,z] = pivot_matlab(B,A,b,C,z,L,e); %换基，即进行初等行变换
        end
        iters = iters + 1;
    end
    x = zeros(n+m,1);                               % [13,16].
    x(B) = b; 
end
