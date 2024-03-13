function [u] = wls_ca_4df_pv_limit(v, u, p_limits, v_limits)
% function [u] = wls_alloc_mch(v, u, umin, umax)
%  [u] = wls_alloc_mch(v,u,p_limits,v_limits)
% WLS_ALLOC - Control allocation using weighted least squares.
%  [u,W,iter] = wls_alloc(B,v,umin,umax,[Wv,Wu,ud,gamma,u0,W0,imax])
% Solves the weighted, bounded least-squares problem
%   min ||Wu(u-ud)||^2 + gamma ||Wv(Bu-v)||^2
%   subj. to  umin <= u <= umax
% using an active set method.
%  Inputs:
%  -------
% v     commanded virtual control (k x 1)
% u0    initial point (m x 1)
% W0    initial working set (m x 1) [empty]
% imax  max no. of iterations [100]
%  Outputs:
%  -------
% u     optimal control
% W     optimal active set
% iter  no. of iterations (= no. of changes in the working set + 1)
%                            0 if u_i not saturated
% Working set syntax: W_i = -1 if u_i = umin_i
%                           +1 if u_i = umax_i
% B         control effectiveness matrix (k x m).
% umin      lower position limits (m x 1).
% umax      upper position limits (m x 1).
% Wv        virtual control weighting matrix (k x k) [I].
% Wu        control weighting matrix (m x m) [I].
% gam       gamma weight (scalar) [1e6].
% ud        desired control (m x 1) [0].
% imax      maximum iterations.
% See also: WLSC_ALLOC, IP_ALLOC, FXP_ALLOC, QP_SIM.
% param of ducted fan
%% 弧度单位
%==============使用等价模型，弧度单位==============
B=[-0.5   0       0.5   0;
    0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
%====仅幅值约束================
if(~v_limits)
umin=[1;1;1;1]*(-p_limits)*pi/180;
umax=[1;1;1;1]*p_limits*pi/180;
else
%====幅值、速度约束================ 
umin=max([1;1;1;1]*(-p_limits)*pi/180,-0.01*400*pi/180+u);
umax=min([1;1;1;1]*p_limits*pi/180,0.01*400*pi/180+u);
end
%%
%========计算当前有效集============
W=zeros(4,1);
infeasible1 = (u <= umin);
W(infeasible1)=-1;
infeasible2 = (u >= umax);
W(infeasible2)= 1;
%===============期望舵机位置========================
ud=[0;0;0;0];% 可作为输入
% Number of variables.
m = 4;
%================================
% u0 = (umin+umax)/2;
% W0 = zeros(m,1);
% 加权系数
gam=1e6;
% 加权矩阵
Wv=eye(3);
Wu=eye(4);
%Wv1=Wv*0.5*y;
% 迭代次数上限
imax=100;
u = wls_alloc_gen(B,v,umin,umax,Wv,Wu,ud,gam,u,W,imax,m);

    
end