clc;clear all;
close all;
addpath(genpath(pwd))
%% setup aircraft and load input data
B=[-0.5     0       0.5     0;
     0      -0.5     0       0.5;
     0.25    0.25    0.25    0.25];
[k,m] = size(B);
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
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
m1=[0.2;0.1;0];
[u, errout, lambda] = DP_LPCA([0;0;0.2],m1,B,umin,umax,100,1)
 m_real=B*u
(m_real-m1)/lambda
if(errout~=0)

[u1, errout1, lambda1] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1)
 m_real1=B*u1
 m_real1/lambda1
end


 %% test for the diffrent lambda behavior on unatainable moment
% [u2, errout2, lambda2] = DP_LPCA([0.5;0;0],[0;0;0],B,umin,umax,100,1)
% [u3, errout3, lambda3] = DP_LPCA([0.5;0;0],[0;0;0],B,umin,umax,100,2000)
% is the same!!!
%% for DP_LPCA
% m可达。lamb=1
% m不可达，0<=lamb<1s, lamb=0对应m=∞

%% for Prioritizing Commands in A.5.1 Dual Branch based DP_LPCA and 0<=lamb<=1
% m1+m2可达。lamb=1
% m1可达,m2不可达，0<=lamb<1, lamb=0对应m2=∞
% 无解，则迭代：删除m2,构造新问题。min-lambda, s.t. Bu=lambda*m1 and 0<=lamb<=1, u_l<=u<=u_u.
% 多个优先级则是上面规则的套娃。

