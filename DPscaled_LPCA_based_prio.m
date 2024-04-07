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
tol = 1e-7;
% then we can set lambda = 1/tol.
%% 
m1=[0.1;0.5;0];
% [u,~,errout,lambda] = DPscaled_LPCA(m1,B,umin,umax,100)
%  m_real=B*u
% (m_real)/lambda
%% and follow method is the same as DPscaled_LPCA
% [u1, errout1, lambda1] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1/tol)
%  m_real1=B*u1
% (m_real1)/lambda1
%% 
m2=[0;0;0.2];
[u2, errout2, lambda2] = DP_LPCA(m2,m1,B,umin,umax,100,1/tol)
if(errout2~=0) % 构造新问题
    [u3, errout3, lambda3] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1/tol)
    m_real3=B*u3
    m_real3/lambda3
else
    m_real2=B*u2
    m_real2*lambda2
    % but the correct is 
    (m_real2*lambda2-m1)/lambda2+m1

end
%% for DPscaled_LPCA
% m可达。lambda对应拉伸到达边界的缩放因子。
% m不可达，lambda对应衰减到达边界的缩放因子。
%% for DP_LPCA and 0<=lambda but not upper bound 
% m可达。lambda对应拉伸到达边界的缩放因子。
% m不可达，lambda对应衰减到达边界的缩放因子。
%%  for Prioritizing Commands  based DP_LPCA but not upper bound
% m1+m2可达。
% m1可达,m2不可达，
% 无解，