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
tol = eps;
% then we can set lambda = 1/tol.
%% 
m1=[0.0;0.0;0.5];
% [u,~,errout,lambda] = DPscaled_LPCA(m1,B,umin,umax,100)
%  m_real=B*u
% (m_real)/lambda
%% and follow method is the same as DPscaled_LPCA
% [u1, errout1, lambda1] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1/tol)
%  m_real1=B*u1
% (m_real1)/lambda1
%% 
m2=[0.1;0.1;0.1];
disp('原问题解：');
[u2, errout2, lambda2] = DP_LPCA(m2,m1,B,umin,umax,100,1/tol)
B*u2
if(errout2~=0) % No Initial Feasible Solution found
    % 构造新问题
    disp('无解，m1+m2无法通过缩放m2穿过边界.隐含m1不可达'); 
    disp('构造新问题'); 
    [u3, errout3, lambda3] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1/tol)
    m_real3=B*u3
    m_real3/lambda3
    disp('或者采取m1+m2的保方向分配，这种情况对应没有更多优先级分解的量'); 
    [u_all, errout_all, lambda_all] = DP_LPCA(m2+m1,[0;0;0],B,umin,umax,100,1/tol)
    disp('最终力矩：');
    B*u_all
else % get a feasible solultion.  % 利用m1可不可达划分
    disp('有解，m1+m2通过缩放m2穿过边界');
    disp('计算m1');
    [u_m1, errout_m1, lambda_m1] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1/tol)
    if(lambda_m1<1)  
        disp('m1不可达,采取m1+m2的保方向分配'); % 包含m1不可达但m1+m2可达的情况，还包括其他仅仅是缩放穿过的情况
        [u_all, errout_all, lambda_all] = DP_LPCA(m2+m1,[0;0;0],B,umin,umax,100,1/tol)
        disp('最终力矩：');
        B*u_all
    else
        disp('m1可达');
        disp('需要修正');
        u=(u2*lambda2-u_m1)/lambda2+u_m1
        B*u
    end 


end
%% for DPscaled_LPCA
% m可达。lambda对应拉伸到达边界的缩放因子。
% m不可达，lambda对应衰减到达边界的缩放因子。
%% for DP_LPCA and 0<=lambda but not upper bound 
% m可达。lambda对应拉伸到达边界的缩放因子。
% m不可达，lambda对应衰减到达边界的缩放因子。
%%  for Prioritizing Commands  based DP_LPCA but not upper bound
% m1+m2可达。修正: 从m1+m2方向收缩。与从m2方向收缩效果一样。
% m1可达,m2不可达（m1+m2不可达），在m2和边缘交点处截断。求解原式得到最佳解。
% m1不可达，修正。此时m2可以拉长，可以得到穿透边缘的解。

% tip：1.多优先级都可以归纳为m1+m2的情形。2.不考虑m1不可达的问题可以使问题简单。
% lambda上界为1时，多种情形在原式下即可表达，不需要修正。但lambda的上界倾向于分配饱和的结果。
% lambda无上界时，需要修正

%% => 由上面分析，逻辑可以是
% 先判断与否，若可达则从m1+m2方向收缩，否则下一步
% m1+m2不可达，m1若可达，则结果自然正确，在m2和边缘交点处截断。
% m1不可达，构造新问题。
%% 以上就是以前的方案。重新整理如下
% [u_all, errout_all, lambda_all] = DP_LPCA(m1+m2,[0;0;0],B,umin,umax,100,1/tol)
% % the problem above always have solution
% if(lambda_all<1)
%     disp('m1+m2不可达, 计算m1可达性');
%     [u_m1, errout_m1, lambda_m1] = DP_LPCA(m1,[0;0;0],B,umin,umax,100,1/tol)
%     if(lambda_m1<1)
%         disp('m1不可达, 重新构造问题');
% 
%     else 
%         disp('m1可达，利用原始问题');
%         [u, errout, lambda] = DP_LPCA(m2,m1,B,umin,umax,100,1/tol)
%     end
% 
% else
%     disp('m1+m2可达, 从m1+m2方向收缩');
%     disp('最终分配结果：');
%     u_all
%     disp('最终力矩：');
%     B*u_all
%     disp('如果从m2方向收缩，会出现不符合要求对值，特别是m1不可达时');
% 
% end


% 不考虑m1不可达的问题可以使问题简单。即m1可达，利用原始问题即为所提出优先级分配


