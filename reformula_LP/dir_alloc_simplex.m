function [u,z,iters] = dir_alloc_simplex(B, v, umin,umax)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=C*x   subj. to  A*x (=、 >=、 <=) b
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
%   x                       [I ;-I ]x <= [umax; Inf; -umin; 0]
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
[k,m]=size(B);

Aeq=[B -v]; % k x (m+1)
beq=zeros(k,1); % k x1
G=[eye(m+1);-eye(m+1)]; % 2(m+1) x (m+1)
h=[umax; 1e4; -umin; 0]; % 20 is the max value of a  % 2(m+1) x1
b=[beq;h];
%% 构造线性规划标准型
% Convert free variables to positively constrained variables
Ad=[Aeq -Aeq; G -G]; % k+2(m+1) x 2(m+1)

% 先把前三个等式的基找到，并化简
P=[Ad(1:k,1:k) zeros(k,2*(m+1));Ad(k+1:k+2*(m+1),1:k) eye(2*(m+1))];
Ad_eye=P\Ad;% 化简 P_inv=inv_mch(P);Ad_eye=P_inv*Ad;% 抖动原因
A=[Ad_eye(1:k,1:2*(m+1)) zeros(k,2*(m+1)); Ad_eye(k+1:k+2*(m+1),1:2*(m+1)) eye(2*(m+1))];
[mad,nad]= size(A);% [mad,~]= size(Ad); %mad=k+2*(m+1)
c =[zeros(1,m) -1]; 
C=[c -c zeros(1,2*(m+1))]; % 1 x 2(m+1)
basis=[1:k 2*(m+1)+1:2*(m+1)+2*(m+1)];% 转C需要特别注意下标的区别 % k
%% Simplex algorithm
%% Iterate through simplex algorithm main loop
z=0;[x,z,iters]=Simplex_loop_matlab(basis, A, b, C, z); % a上界小时会抖动 % xx_matlab come from xx_mch
% [x,z,~,~]=Simplex_loop_C(basis, A, b, C,mad,nad); % 线性规划单纯形法
% 转化解
u1=x(1:m)-x(m+1+1:m+1+m);
if z>1  % 放大了倍数，再还原，若小于1，则表示需要缩小，x已经自然到达边界
    u = u1./(z);
else
    u=u1;
end
end