clear all;
close all;
addpath(genpath(pwd))

%% setup aircraft and load input data
% B=[-0.5     0       0.5     0;
%      0      -0.5     0       0.5;
%      0.25    0.25    0.25    0.25];
l1=0.148;l2=0.069;k_v=3;
B=k_v*[-l1     0       l1     0;
     0      -l1     0       l1;
     l2    l2    l2    l2];
[k,m] = size(B);
% u_0=ones(m,1)*20*pi/180;
u_0=[0.0122;
     0.0122;
     0.3491;
     0.3491];
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
% plim=[umin-u_0 umax-u_0];
% q=vview(B,plim,pinv(B));
% run Generate_input_data;
load 'input.mat'; % get v and the len_command_px4 (len_command_px4 is size of command_px4, which come from flght log data)
[~,N]=size(v);
%% setup function of allocation lib
% ========
%% setup ACA
% global NumU
NumU=m;
LPmethod=3; % LPmethod should be an integer between 0 and 5. when LPmethod=2 set upper of lambda to Inf can't save this method!!! but big number is the same as that based linprog
INDX=ones(1,m);  % active effectors
IN_MAT = [B     zeros(k,1)
          umin' 0
          umax' 0
          INDX  LPmethod];
IN_MAT1 = [B     zeros(k,1)
          (umin-u_0)' 0
          (umax-u_0)' 0
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
%%
u=zeros(m,1);
x_LPwrap=zeros(m,N);
x_LPwrap_incre=zeros(m,N);
x_allocator_dir_LPwrap_4=zeros(m,N);
x_CGIwrapp=zeros(m,N);
x_DAwrap=zeros(m,N);
x_VJAwrap=zeros(m,N);
x_inv=zeros(m,N);
x_wls=zeros(m,N);
x_wls_gen=zeros(m,N);
x_dir_alloc_linprog=zeros(m,N);
x_dir_alloc_linprog_re=zeros(m,N);
x_dir_alloc_linprog_re_bound=zeros(m,N);
x_use_LP_lib=zeros(m,N);
tic;
%% simulate flight process  
for idx=1:N
    
    % IN_MAT(1:3,end) = v(:,idx);

    u = LPwrap(IN_MAT,v(:,idx),NumU); % function of ACA lib
    x_LPwrap(:,idx) = min(max(u, umin), umax);

    % u = LPwrap(IN_MAT1,v(:,idx),NumU); % incremental form. 0 have to be a feasible solution.
    % x_LPwrap_incre(:,idx) = min(max(u, umin-u_0), umax-u_0)+u_0;

    % u= CGIwrap(IN_MAT,v(:,idx),NumU);
    % x_CGIwrapp(:,idx) = min(max(u, umin), umax);

    % u = DAwrap(IN_MAT,v(:,idx),NumU);
    % x_DAwrap(:,idx) = min(max(u, umin), umax);

    % u = VJAwrap(IN_MAT,v(:,idx),NumU);
    % x_VJAwrap(:,idx) = min(max(u, umin), umax);

    % u=pinv(B)*v(:,idx);
    % x_inv(:,idx) = min(max(u, umin), umax);
    % 
    % [u,~,~] = wls_alloc(B,v(:,idx),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
    % x_wls(:,idx) = min(max(u, umin), umax);
    % 
    % u =wls_alloc_gen(B,v(:,idx),umin,umax,eye(k),eye(m),zeros(m,1),1e6,zeros(m,1),zeros(m,1),100,4);
    % x_wls_gen(:,idx) = min(max(u, umin), umax);
    % 
    % [u,~] = dir_alloc_linprog(B,v(:,idx), umin, umax, 1e4); % LPmethod=2 and lam=1 of dir_alloc_linprog is lager but similar
    % x_dir_alloc_linprog(:,idx) = min(max(u, umin), umax);
    % 
    % [u,~] = dir_alloc_linprog_re(B,v(:,idx), umin, umax);
    % x_dir_alloc_linprog_re(:,idx) = min(max(u, umin), umax);
    % 
    % [u,~] = dir_alloc_linprog_re_bound(B,v(:,idx), umin, umax, 1e4);% the
    % same as dir_alloc_linprog for any lam >=1, lam have to be >1 when use
    % linprog, that will be the same as LPmethod=3
    % x_dir_alloc_linprog_re_bound(:,idx) = min(max(u, umin), umax);
    % 
    % [u,~] = use_LP_lib(B,v(:,idx), umin, umax); % ToDo: use the LP lib
    % x_use_LP_lib(:,idx)=min(max(u, umin), umax);
    % 
    % [u,~,~] =allocator_dir_LPwrap_4(single(B), single( v(:,idx)), single(umin),single(umax)); % ToDo: 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.914283e-09。 
    % x_allocator_dir_LPwrap_4(:,idx) = min(max(u, umin), umax);
    
end
elapsed_time = toc;
fprintf('代码执行时间：%.2f 秒\n', elapsed_time);
%% Determine the variables to use for comparison.
% run target of alloc_cpp (./main) to generate output.csv
output = readmatrix('output.csv')';% or delete this line to just compare the matlab implement method
command_px4=v(:,1:len_command_px4);
% just use the flight data to compare.

x1=output(:,1:len_command_px4); % or x_xxx above
% x1=x_LPwrap_incre(:,1:len_command_px4);
x2=x_LPwrap(:,1:len_command_px4);

% actual moments produced. The B matrix have to be the same.
U1=B*x1;
U2=B*x2;

dt=mean(delta_t_s);
t=0:dt:dt*(len_command_px4-1);
tt=1:1:(len_command_px4-1);

error1=U1-command_px4;
error2=U2-command_px4;
figure,
subplot(4,1,1)
plot(t,x1(1,:),'r-');hold on;
plot(t,x2(1,:),'b--');hold on;
% plot(t,u_px4(:,1),'g.');hold on;
subplot(4,1,2)
plot(t,x1(2,:),'r-');hold on;
plot(t,x2(2,:),'b--');hold on;
% plot(t,u_px4(:,2),'g.');hold on;
subplot(4,1,3)
plot(t,x1(3,:),'r-');hold on;
plot(t,x2(3,:),'b--');hold on;
% plot(t,u_px4(:,3),'g.');hold on;
subplot(4,1,4)
plot(t,x1(4,:),'r-');hold on;
plot(t,x2(4,:),'b--');hold on;
% plot(t,u_px4(:,4),'g.');hold on;
% figure,
% subplot(3,1,1)
% plot(t,error1(1,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
% plot(t,error2(1,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % plot(t,error1(1,:)-error2(1,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
% % plot(t,U1(1,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
% % plot(t,U2(1,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % plot(t,command_px4(:,1),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
% 
% subplot(3,1,2)
% plot(t,error1(2,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
% plot(t,error2(2,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % plot(t,error1(2,:)-error2(2,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
% % plot(t,U1(2,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
% % plot(t,U2(2,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % plot(t,command_px4(:,2),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
% subplot(3,1,3)
% plot(t,error1(3,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
% plot(t,error2(3,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % plot(t,error1(3,:)-error2(3,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
% % plot(t,U1(3,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
% % plot(t,U2(3,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % plot(t,command_px4(:,3),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;

% outside_x1=output(:,len_command_px4+1:end);
% outside_x1=x_LPwrap_incre(:,len_command_px4+1:end);
% outside_x2=x_LPwrap(:,len_command_px4+1:end);
% outside_U1=B*outside_x1;
% outside_U2=B*outside_x2;
% outside_err=outside_x1-outside_x2;
% figure,
% plot3(outside_U1(1,:),outside_U1(2,:),outside_U1(3,:),'r*');hold on;
% % figure,
% plot3(outside_U2(1,:),outside_U2(2,:),outside_U2(3,:),'g*');

