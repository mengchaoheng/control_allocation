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
m1=[0.1;0.1;0.0]; % [0;0;0.2]; or [0;0;0.5]; % higher
m2=[0.0;0.0;0.0];% [0.1;0.1;-0.4]; or [0.1;0.1;0.4]; % lower
disp('原问题解：');
% tic;
[u, errout, lambda] = DP_LPCA_copy1(m1,m2,B,umin,umax,100)
% toc;
u=restoring_cpp(B,u,umin,umax)
B*u
add=m1+m2
% not work


% when m1_higher is not attainable, No Initial Feasible Solution found,
% resukt is error.

% Or when m2=0, Solution is unbounded, Solver error.