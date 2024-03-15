%% Iterate through simplex algorithm main loop
function [x,z,iters,f]=Simplex_loop_C(B, A, b, c,m,n)
% (c) mengchaoheng
% 不考虑无解的情形
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=、 >=、 <=) b
%   x 
    %% Initialization
    tol=single(1e6);   
    x = single(zeros(n,1));
    iters=uint32(0);
    z=single(0);
    f=false;
%     [m,n] = size(A);
%     while ~all(c>=0)                      % 3.~isempty(c(c(N)<0))
%     e = find(c < 0, 1, 'first'); % 进基变量索引    % 4. e = N(find(c(N)<0,1))
    while(1)
        flag=false;
        e=uint8(0);
        L=uint8(0);
        for i=1:n
            if c(i)< -1/tol % <0
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
            delta=tol*single(ones(m,1));
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
                [B,A,b,c,z] = pivot_C(B,A,b,c,z,L,e,m,n); %换基，即进行初等行变换
            end
            iters = iters + uint32(1);
        else
            f=true;
            break;% 已经找到解，或者无解
        end
    end
    x(B) = b; 
end
