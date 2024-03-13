%% Iterate through simplex algorithm main loop
function [x,z,iters,f]=Simplex_loop_C(B, A, b, C, z,m,n)
% (c) mengchaoheng
% 不考虑无解的情形
% Last edited 2019-11
%   min z=C*x   subj. to  A*x (=、 >=、 <=) b
%   x 
    %% Initialization
    tol=1e6;   
    x = zeros(n,1);
    iters=0;
    z=0;
    f=0;
%     [m,n] = size(A);
%     while ~all(C>=0)                      % 3.~isempty(C(C(N)<0))
%     e = find(C < 0, 1, 'first'); % 进基变量索引    % 4. e = N(find(C(N)<0,1))
    while(1)
        flag=false;
        e=0;
        L=0;
        for i=1:n
            if C(i)< -1/tol % <0
                flag=true;
                e=i;
                break;
            end
        end
        if flag==true
%             a_ie=A(:,e);
%             ip=a_ie>(1/tol);
%             delta=tol*ones(m,1);
%             if ~isempty(ip)
%                 delta(ip)=b(ip)./a_ie(ip);
%             end
            MIN=tol;
            delta=tol*ones(m,1);
            for i=1:m
                a_ie=A(i,e);
                if a_ie>(1/tol)
                    delta(i)=b(i)./a_ie;
                else
                    delta(i)=tol;
                end
                if delta(i)<MIN
                    L=i;
                    MIN=delta(i);
                end
            end
%             [~,L]=min(delta);%选择离基 (离基在B数组中的行索引)
    %         li = B(L);    % 离基变量索引               
%             if delta(L) >= tol    
            if MIN >= tol
                return;  % 无界解
            else % 此时一定有一个L
                [B,A,b,C,z] = pivot_C(B,A,b,C,z,L,e,m,n); %换基，即进行初等行变换
            end
            iters = iters + 1;
        else
            f=1;
            break;% 已经找到解，或者无解
        end
    end
    x(B) = b; 
end
