%% Iterate through simplex algorithm main loop
function [x,z,iters,f]=Simplex_loop_C(B, A, b, c, z,m,n)
% (c) mengchaoheng
% �������޽������
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=�� >=�� <=) b
%   x 
    %% Initialization
    tol=1e6;   
    x = zeros(n,1);
    iters=0;
    z=0;
    f=0;
%     [m,n] = size(A);
%     while ~all(c>=0)                      % 3.~isempty(c(c(N)<0))
%     e = find(c < 0, 1, 'first'); % ������������    % 4. e = N(find(c(N)<0,1))
    while(1)
        flag=false;
        e=0;
        L=0;
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
%             [~,L]=min(delta);%ѡ����� (�����B�����е�������)
    %         li = B(L);    % �����������               
%             if delta(L) >= tol    
            if MIN >= tol
                return;  % �޽��
            else % ��ʱһ����һ��L
                [B,A,b,c,z] = pivot_C(B,A,b,c,z,L,e,m,n); %�����������г����б任
            end
            iters = iters + 1;
        else
            f=1;
            break;% �Ѿ��ҵ��⣬�����޽�
        end
    end
    x(B) = b; 
end
