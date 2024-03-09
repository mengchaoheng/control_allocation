function [u] = dir_alloc_mch(v, umin,umax)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=�� >=�� <=) b
%   x 
% ԭ����
% Performs direct control allocation by solving the LP
%   max z=a   subj. to  Bu = av
%   a,u               umin <= u <= umax
% If a > 1, set u = u/a.
% Note: This function has not been optimized for speed.
%  Inputs:
%  -------
% B     control effectiveness matrix (k x m)
% v     commanded virtual control (k x 1)
% umin  lower position limits (m x 1)
% umax  upper position limits (m x 1)
%  Outputs:
%  -------
% u     optimal control (m x 1)
% a     scaling factor  
%% �����
%   min z=[0 -1]x   subj. to  [B -v]x = 0
%   x                       [I 0;-I 0]x <= [umax; -umin]
%   ���� x=[u; a]
% ��Ӧ��͹�Ż���p139,��Ϊ
%   min z=c*x   subj. to  Aeq*x = beq
%   x                     G*x <= h
% �ϲ�
%   min z=c*x   subj. to  [Aeq; G]*x (=��<=) [beq;h]
%   x                     
% ��֤x>=0������
%   min z=[c -c]*X   subj. to  [Aeq -Aeq;G -G]*X (=��<=) [beq;h]
%    X                                          
% ���� X=[x^+; x^-]
%%
% B=[-0.5   0       0.5   0;
%      0  -0.5    0       0.5;
%     0.25   0.25   0.25   0.25];
% B=[-0.5   0       0.5   0;
%      0  -0.5    0       0.5;
%     0.25   0.25   0.25   0.25];
% Aeq=[B -v];
% beq=zeros(3,1);
% G=[eye(5);-eye(5)];
% h=[umax; 20; -umin; 0];
% ������Թ滮
% b=[beq;h];
b=[zeros(3,1);umax; 20; -umin; 0];
%% �������Թ滮��׼��
% Convert free variables to positively constrained variables
% Ad=[Aeq -Aeq; G -G];
% Ad=[B -v -B v; eye(5) -eye(5);-eye(5) eye(5)];
% Ad=[-0.5000         0    0.5000         0         0    0.5000         0   -0.5000         0         0;
%           0   -0.5000         0    0.5000         0         0    0.5000         0   -0.5000         0;
%      0.2500    0.2500    0.2500    0.2500         0   -0.2500   -0.2500   -0.2500   -0.2500   	    0;
%      1.0000         0         0         0         0   -1.0000         0         0         0         0;
%           0    1.0000         0         0         0         0   -1.0000         0         0         0;
%           0         0    1.0000         0         0         0         0   -1.0000         0         0;
%           0         0         0    1.0000         0         0         0         0   -1.0000         0;
%           0         0         0         0    1.0000         0         0         0         0   -1.0000;
%     -1.0000         0         0         0         0    1.0000         0         0         0         0;
%           0   -1.0000         0         0         0         0    1.0000         0         0         0;
%           0         0   -1.0000         0         0         0         0    1.0000         0         0;
%           0         0         0   -1.0000         0         0         0         0    1.0000         0;
%           0         0         0         0   -1.0000         0         0         0         0    1.0000];
% Adֻ�е�5����10�и���v��ͬ����ͬ�������̶�����
Ad5=[0 0 0 0 0 0 0 1 0 0 0 0 -1]';  
Ad10=[0 0 0 0 0 0 0 -1 0 0 0 0 1]';  
Ad5(1:3) =-v;
Ad10(1:3) =v;
% [mad,~]= size(Ad);
% mad=13;
% �Ȱ�ǰ������ʽ�Ļ��ҵ���������
% B_inv=[Ad(1:3,1:3) zeros(3,mad-3);Ad(4:mad,1:3) eye(mad-3)];
% B_inv=[Ad(1:3,1:3) zeros(3,10);Ad(4:mad,1:3) eye(10)];
% P=[Ad(1:3,1:3) zeros(3,10);Ad(4:13,1:3) eye(10)];
% ����
% Ad_eye=P\Ad;% ����
% �޹��е�����P_inv�ǳ�����
P_inv=[-1     1     2     0     0     0     0     0     0     0     0     0     0;
     0    -2     0     0     0     0     0     0     0     0     0     0     0;
     1     1     2     0     0     0     0     0     0     0     0     0     0;
     1    -1    -2     1     0     0     0     0     0     0     0     0     0;
     0     2     0     0     1     0     0     0     0     0     0     0     0;
    -1    -1    -2     0     0     1     0     0     0     0     0     0     0;
     0     0     0     0     0     0     1     0     0     0     0     0     0;
     0     0     0     0     0     0     0     1     0     0     0     0     0;
    -1     1     2     0     0     0     0     0     1     0     0     0     0;
     0    -2     0     0     0     0     0     0     0     1     0     0     0;
     1     1     2     0     0     0     0     0     0     0     1     0     0;
     0     0     0     0     0     0     0     0     0     0     0     1     0;
     0     0     0     0     0     0     0     0     0     0     0     0     1];
% P_inv=inv_mvh(P);
% Ad_eye=P_inv*Ad;
% Ad_eye=[1     0     0     1     0    -1     0     0    -1     0;
%         0     1     0    -1     0     0    -1     0     1     0;
%         0     0     1     1     0     0     0    -1    -1     0;
%         0     0     0    -1     0     0     0     0     1     0;
%         0     0     0     1     0     0     0     0    -1     0;
%         0     0     0    -1     0     0     0     0     1     0;
%         0     0     0     1     0     0     0     0    -1     0;
%         0     0     0     0     0     0     0     0     0     0;
%         0     0     0     1     0     0     0     0    -1     0;
%         0     0     0    -1     0     0     0     0     1     0;
%         0     0     0     1     0     0     0     0    -1     0;
%         0     0     0    -1     0     0     0     0     1     0;
%         0     0     0     0     0     0     0     0     0     0];
% �������Ϸ�����P_inv*Adֻ�е�5����10����v��ͬ���仯��
% for i=1:13
%     temp1=0;
%     temp2=0;
%     for k=1:13
%         temp1=temp1 + P_inv(i,k)*Ad(k,5);
%         temp2=temp2 + P_inv(i,k)*Ad(k,10);
%     end
%     Ad_eye(i,5)=temp1;
%     Ad_eye(i,10)=temp2;
% end
% �����ɳڱ�����Ӧ�Ļ�
% A=[Ad_eye(1:3,1:10) zeros(3,10); Ad_eye(4:13,1:10) eye(10)];
% A��Ad_eye�����䣬��5����10����P_inv*Ad�仯�Ĳ������йأ������ǳ���
A=[1     0     0     1     0    -1     0     0    -1     0     0     0     0     0     0     0     0     0     0     0;
   0     1     0    -1     0     0    -1     0     1     0     0     0     0     0     0     0     0     0     0     0;
   0     0     1     1     0     0     0    -1    -1     0     0     0     0     0     0     0     0     0     0     0;
   0     0     0    -1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0;
   0     0     0     1     0     0     0     0    -1     0     0     1     0     0     0     0     0     0     0     0;
   0     0     0    -1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0;
   0     0     0     1     0     0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0;
   0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0;
   0     0     0     1     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0;
   0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0;
   0     0     0     1     0     0     0     0    -1     0     0     0     0     0     0     0     0     1     0     0;
   0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0;
   0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1];
for i=1:13
    temp1=0;
    temp2=0;
    for k=1:13
        temp1=temp1 + P_inv(i,k)*Ad5(k);
        temp2=temp2 + P_inv(i,k)*Ad10(k);
    end
    A(i,5)=temp1;
    A(i,10)=temp2;
end

c=[0 0 0 0 -1 0 0 0 0 1 zeros(1,10)];
basis=[1:3 11:20];% תC��Ҫ�ر�ע���±������
z = 0;
%% Simplex algorithm
%% Iterate through simplex algorithm main loop
% [x,z]=Simplex_loop_mch(basis, A, b, c, z); % ���Թ滮�����η�
%% Initialization
iters=0;
e=0;
% [m,n] = size(A);
% m=13;
% n=20;
while ~all(c>=0)                      % 3.~isempty(c(c(N)<0))
%     e = find(c < 0, 1, 'first'); % ������������    % 4. e = N(find(c(N)<0,1))
    for i=1:20
        if(c(i)<0)
            e=i;
            break;
        end
    end
    a_ie=A(:,e);
    ip=a_ie>0;
    delta=1e6*ones(13,1);
    if ~isempty(ip)
        delta(ip)=b(ip)./(a_ie(ip));
    end
    [~,L]=min(delta);%ѡ����� ������
%         li = basis(L);    % �����������               
    if delta(L) == 1e6                         
        % disp('System is unbounded!');
        break;  % 11.
    else
        %% Perform pivot operation, exchanging L-row with e-coLumn variabLe
%         [basis,A,b,c,z] = pivot_mch(basis,A,b,c,z,L,e);
%         %�����������г����б任��תCҪ�ر�ע�⣬Ҫ��ϵ���ȱ����ٽ���forѭ��
        % Compute the coefficients of the equation for new basic variabLe x_e.
        a_Le=A(L,e);% 4.    
        b(L) = b(L)/a_Le;                 
            
            A(L,1:20) = A(L,1:20)./(a_Le);            
            % Compute the coefficients of the remaining constraints.
            i=[1:L-1 L+1:13];     %  i = find(B~=li);
            if ~isempty(i)
                b(i) = b(i) - A(i,e)*b(L);                                               
                A(i,1:20) = A(i,1:20) - A(i,e)*A(L,1:20);	
            end
            % Compute the objective function
            z = z - c(1,e)*b(L);                     
            c(1,1:20) = c(1,1:20) - c(1,e) * A(L,1:20);          
            % Compute new sets of basic and nonbasic variabLes.
            % N(find(N==e,1)) = li; 
            basis(L)=e;                 %  B(find(B==li,1)) = e;
    end
    iters = iters + 1;
end
x = zeros(20,1);                               % [13,16].
x(basis) = b; 
% ת����
u1=x(1:4)-x(6:9);
if z>1  % �Ŵ��˱������ٻ�ԭ����С��1�����ʾ��Ҫ��С��x�Ѿ���Ȼ����߽�
    u = u1./(z);
else
    u=u1;
end
end