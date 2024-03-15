clc;
clear all;
close all;
addpath(genpath(pwd))
folder ='some_modified_function'; 
rmpath(folder) % remove old version
folder ='s-function_used_in PlanD'; 
rmpath(folder) % remove old version
%%

% Simulation:
%   qp_sim        - Response of static QP allocator.
%   dir_sim       - Response of direct allocator.
%   dyn_sim       - Response of dynamic allocator.
% Then call generic control allocation simulation subroutine.
%   alloc_sim     - Control allocation simulation.
% Then call:
% QP based control allocation
%     sls_alloc     - Active set, sequential least squares.
%     wls_alloc     - Active set, weighted least squares.
%     wlsc_alloc    - C implementation of wls_alloc.
%     mls_alloc     - Active set, minimal least squares.
%     ip_alloc      - Interior point method.
%     cgi_alloc     - Cascaded generalized inverses method.
%     fxp_alloc     - Fixed-point method.
% 
% Direct allocation.
%     dir_alloc     - Direct control allocation.


% ========
% test by the lib code
% F18 example
load f18data
% qp_sim
% [u,W,time,iter] = qp_sim(B,v,plim,varargin) 
% [u,W,time,iter] = qp_sim(B,v,plim,[rlim,T,Wv,Wu,ud],options)
%  Options: options = option_1,value_1,option_2,value_2,...
%  --------
% 'alg'    numerical algorithm: 'sls'    SLS_ALLOC
%                               'wls'    WLS_ALLOC (default)
%                               'wlsc'   WLSC_ALLOC
%                               'mls'    MLS_ALLOC
%                               'ip'     IP_ALLOC
%                               'cgi'    CGI_ALLOC
%                               'fxp'    FXP_ALLOC
% 'imax'   max no. of iterations [100]
% 'gam'    weight used in algorithms based on weighted LS [1e6]
% 'tol'    tolerance used in IP_ALLOC stopping criterion [1e-6]
% 'hot'    hotstart solver (not ip/cgi) with previous solution (0/[1])
% 'ui'     initial control signal
% 'Wi'     initial active constraints
% 'rep'    no. of repetitions [1]
[u,W,time,iter]=qp_sim(B,v,plim,rlim,T1,'alg','sls'); % change alg for other method test
figure,
plot(tn,u);
figure,
plot(tn,v,'k',tn,B*u);
% dir_sim
u=dir_sim(B,v,plim,rlim,T2);
figure,
plot(tn,u*180/pi),ylabel('Controls (deg)')
figure,
plot(tn,B*u,tn,v,'k--'),legend('roll','pitch','yaw')
%  [u,W,time,iter] = dyn_sim(B,v,plim,[rlim,T,Wv,W1,W2,S],options)
%  Step response example:
B = [2 1]; t = 0:.2:10; v = 1*(t>1); plim = [-1 1;-1 1];
W1 = eye(2); W2 = diag([5 0]); S = pinv(B);
u = dyn_sim(B,v,plim,[],.2,1,W1,W2,S);
figure,
bodemag(dca(B,S,W1,W2,.2))
figure,
stairs(t,[u' v']),legend('u_1','u_2','v=2u_1+u_2')