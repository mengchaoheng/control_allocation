function [u,z,iters] = dir_alloc_sim_new(v, umin,umax,B)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=、 >=、 <=) b
%   x 
% 原问题
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
%% 整理成
%   min z=[0 -1]x   subj. to  [B -v]x = 0
%   x                       [I 0;-I 0]x <= [umax; -umin]
%   其中 x=[u; a]
% 对应《凸优化》p139,记为
%   min z=c*x   subj. to  Aeq*x = beq
%   x                     G*x <= h
% 合并
%   min z=c*x   subj. to  [Aeq; G]*x (=、<=) [beq;h]
%   x                     
% 保证x>=0，变形
%   min z=[c -c]*X   subj. to  [Aeq -Aeq;G -G]*X (=、<=) [beq;h]
%    X                                          
% 其中 X=[x^+; x^-]
%%
% B=[-0.5   0       0.5   0;
%      0  -0.5    0       0.5;
%     0.25   0.25   0.25   0.25];

Aeq=[B -v];
beq=zeros(3,1);
G=[eye(5);-eye(5)];
h=[umax; 20; -umin; 0];
%b求解线性规划
b=[beq;h];

%% 构造线性规划标准型
% Convert free variables to positively constrained variables
Ad=[Aeq -Aeq; G -G];
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
% Ad只有第5，第10列根据v不同而不同，其他固定不变
% [mad,~]= size(Ad);
mad=13;
nad=20;
% 先把前三个等式的基找到，并化简
P=[Ad(1:3,1:3) zeros(3,mad-3);Ad(4:mad,1:3) eye(mad-3)];
% P=[Ad(1:3,1:3) zeros(3,10);Ad(4:mad,1:3) eye(10)];
% P=[Ad(1:3,1:3) zeros(3,10);Ad(4:13,1:3) eye(10)];
% 求逆
% Ad_eye=P\Ad;% 化简
% 无关列的逆阵P_inv是常矩阵
% P_inv=[-1     1     2     0     0     0     0     0     0     0     0     0     0;
%      0    -2     0     0     0     0     0     0     0     0     0     0     0;
%      1     1     2     0     0     0     0     0     0     0     0     0     0;
%      1    -1    -2     1     0     0     0     0     0     0     0     0     0;
%      0     2     0     0     1     0     0     0     0     0     0     0     0;
%     -1    -1    -2     0     0     1     0     0     0     0     0     0     0;
%      0     0     0     0     0     0     1     0     0     0     0     0     0;
%      0     0     0     0     0     0     0     1     0     0     0     0     0;
%     -1     1     2     0     0     0     0     0     1     0     0     0     0;
%      0    -2     0     0     0     0     0     0     0     1     0     0     0;
%      1     1     2     0     0     0     0     0     0     0     1     0     0;
%      0     0     0     0     0     0     0     0     0     0     0     1     0;
%      0     0     0     0     0     0     0     0     0     0     0     0     1];
P_inv=inv_mch(P);
Ad_eye=P_inv*Ad;
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
% 根据以上分析，P_inv*Ad只有第5，第10列依v不同而变化。
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
% 加上松弛变量对应的基
A=[Ad_eye(1:3,1:10) zeros(3,10); Ad_eye(4:13,1:10) eye(10)];
% A是Ad_eye的扩充，第5，第10列与P_inv*Ad变化的部分列有关，其他是常数
% A=[1     0     0     1     0    -1     0     0    -1     0     0     0     0     0     0     0     0     0     0     0;
%    0     1     0    -1     0     0    -1     0     1     0     0     0     0     0     0     0     0     0     0     0;
%    0     0     1     1     0     0     0    -1    -1     0     0     0     0     0     0     0     0     0     0     0;
%    0     0     0    -1     0     0     0     0     1     0     1     0     0     0     0     0     0     0     0     0;
%    0     0     0     1     0     0     0     0    -1     0     0     1     0     0     0     0     0     0     0     0;
%    0     0     0    -1     0     0     0     0     1     0     0     0     1     0     0     0     0     0     0     0;
%    0     0     0     1     0     0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0;
%    0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0;
%    0     0     0     1     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0;
%    0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0     0     1     0     0     0;
%    0     0     0     1     0     0     0     0    -1     0     0     0     0     0     0     0     0     1     0     0;
%    0     0     0    -1     0     0     0     0     1     0     0     0     0     0     0     0     0     0     1     0;
%    0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1];
% for i=1:13
%     temp1=0;
%     temp2=0;
%     for k=1:13
%         temp1=temp1 + P_inv(i,k)*Ad5(k);
%         temp2=temp2 + P_inv(i,k)*Ad10(k);
%     end
%     A(i,5)=temp1;
%     A(i,10)=temp2;
% end

c=[0 0 0 0 -1 0 0 0 0 1 zeros(1,10)];
basis=[1:3 11:20];% 转C需要特别注意下标的区别
z = 0;
%% Simplex algorithm
%% Iterate through simplex algorithm main loop
[x,z,iters]=Simplex_loop_C(basis, A, b, c, z,mad,nad); % 线性规划单纯形法
u1=x(1:4)-x(6:9);
if z>1  % 放大了倍数，再还原，若小于1，则表示需要缩小，x已经自然到达边界
    u = u1./(z);
else
    u=u1;
end
end