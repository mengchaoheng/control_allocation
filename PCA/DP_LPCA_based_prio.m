clc;clear all;
close all;
addpath(genpath(pwd))
%% setup aircraft and load input data
%=============================4==================================
B=[-0.5     0       0.5     0;
     0      -0.5     0       0.5;
     0.25    0.25    0.25    0.25];
%=============================6==================================
% d=60*pi/180;
% B_inv=[-1 0 1;-1 -1 1;1 -1 1;1 0 1;1 1 1;-1 1 1];B=pinv(B_inv);
[k,m] = size(B);
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
plim=[umin umax];
% q=vview(B,plim,pinv(B))
%% setup function of allocation lib
% ========
%% setup ACA
global NumU

NumU=m;
LPmethod=2; % LPmethod should be an integer between 0 and 5. when LPmethod=2 set upper of lambda to Inf can't save this method!!! but big number is the same as that based linprog
INDX=ones(1,m);  % active effectors
IN_MAT = [B     zeros(k,1)
          umin' 0
          umax' 0
          INDX  LPmethod];
%% setup qcat. just wls_alloc and not a Hotstart setting here, use test_qcat.m for more test, 
Wv   = eye(k);     % QP allocation
Wu   = eye(m);
ud   = zeros(m,1);
gam  = 1e6;	     % weight
u0 = (umin+umax)/2;
W0 = zeros(m,1);
imax = 100;	     % no of iterations
% ========
% tol = 1e-7;
% then we can set lambda = 1/tol.
%% simulate flight process   
m1=[0.2;0.2;0.0]; % [0;0;0.2]; or [0;0;0.5]; % higher
m2=[-0;0.0;0.34];% [0.1;0.1;-0.4]; or [0.1;0.1;0.4]; % lower

u_inv=min(max(pinv(B)*(m1+m2), umin), umax)
B*u_inv
disp('原问题解：');
% tic;
[u, errout, lambda] = DP_LPCA_copy(m1,m2,B,umin,umax,100);
% toc;
if(errout<0)
    % 构造新问题
    disp('errout==-2无解，即m1+m2不可达且m1不可达');% 不可通过单独收缩m2实现.此时运行1.5个DP_LPCA_copy，因为第一个只运行了一半，到No Initial Feasible Solution found 即 m_higher unattainable
    disp('or errout==-1 m2=0');% m2=0是直接退出第一个DP_LPCA_copy，此时运行一个DP_LPCA_copy
    [u1, errout1, lambda1] = DP_LPCA_copy([0;0;0],m1+m2,B,umin,umax,100) % second pram is m1 or  m1+m2
    u1=restoring_cpp(B,u1,umin,umax)
    m_real1=B*u1
    (m1+m2)*lambda1
else
    disp('有解，m1可达或者m1+m2可达'); %运行一个DP_LPCA_copy的时间
    u=restoring_cpp(B,u,umin,umax) % 
    disp('原问题解产生的力矩总是满足要求');
    m_real=B*u
    disp('解统一表达式是m1+lambda*m2');
    m1+lambda*m2
    % if(lambda<1)
    %     disp('m1可达（无论m1+m2取何值总能有解），但lamb<1表示m1+m2不可达');
    %     disp('解是m1+lambda*m2:');
    %     m1+lambda*m2
    % else
    %     disp('m1+m2可达（无论m1），lamb=1');
    %     disp('解是m1+m2:');
    %     m1+m2
    % end
end
% disp('DP_LPCA_prio：');
% [u, ~, ~] = DP_LPCA_prio1(m1,m2,B,umin,umax,100)
% u=restoring_cpp(B,u,umin,umax)
%% for DP_LPCA_copy
% m可达。lamb=1
% m不可达，0<=lamb<1s, lamb=0对应m=∞

%% for Prioritizing Commands in A.5.1 Dual Branch based DP_LPCA_copy and 0<=lamb<=1
% m1+m2可达。lamb=1
% m1可达,m2不可达，0<=lamb<1, lamb=0对应m2=∞
% 无解，没有保方向的特性。则迭代：删除m2,构造新问题。min-lambda, s.t. Bu=lambda*m1 and 0<=lamb<=1,
% u_l<=u<=u_u.只要m1不可达一定无解
% 多个优先级则是上面规则的套娃。

