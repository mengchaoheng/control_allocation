


clear all;
close all;clc
addpath(genpath('/Users/mch/Proj/control_allocation'));
%% 全局图形设置
set(groot, ...
    'defaultAxesFontSize', 8, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1);
%% setup aircraft and load input data
mass=0.75;g=9.8;
I=diag([2.5,2.1,4.3]);
L=0.14;
beta=56*pi/180;
c_q=2.37e-8;
c_t=1.51e-6;
% B= [1 1 1 1;L*sin(beta) -L*sin(beta) -L*sin(beta) L*sin(beta);-L*cos(beta) -L*cos(beta) L*cos(beta) L*cos(beta);c_q/c_t -c_q/c_t c_q/c_t -c_q/c_t];
% B= [1  0;0 1]; % for zero
% B= [1  1 0;0 0  1];% for robustly
% B= [1  0 -0.5;0 1  -0.5]; % for case
% B= [1  1 1 1;1 -1 -1 1;-1 -1 1 1;1 -1 1 -1];% for full quad 
B= [1 -1 -1 1;-1 -1 1 1;1 -1 1 -1]; % for tau of quad
% B= [1  1 1 1;1 -1 -1 1;-1 -1 1 1]; % for T and tilt of quad

[k,m] = size(B);
%检查robust
% 枚举所有从 m 列中选出 k 列的组合
comb = nchoosek(1:m, k);    % size(comb) = [C(m,k) × k]

% 预分配结果向量
numComb = size(comb, 1);
rks = zeros(numComb,1);

% 对每个组合计算秩
for i = 1:numComb
    subB = B(:, comb(i,:));   % 取出第 k 种列组合形成的 n×n 子矩阵
    rks(i) = rank(subB);
end

% 如果你只关心哪些子矩阵是满秩（秩 = k），可以这样：
fullRankIdx = find(rks == k);      % 返回所有满秩子矩阵在 comb 中的行号
fullRankCombs = comb(fullRankIdx, :);  % 对应的列索引组合

% 输出一下结果
if numel(fullRankIdx)==numComb
  fprintf('共有 %d 个 %d×%d 子矩阵是满秩的, robust。\n', numel(fullRankIdx), k, k);
else
    fprintf('不满足 robust。\n');
end


umin=ones(m,1)*0.0;
umax=ones(m,1)*8;
plim=[umin umax];

% run Generate_input_data;
% ---- User-level design parameters ----
r0   = 0.00;      % start radius
r1   =  12;      % end   radius (e.g., your r_d)
z0   = 0.00;      % start height
z1   = 7;      % end   height
n    = 5.0;       % desired turns
T    = 10;       % total time duration [s]
dt   = 0.01;     % time step

% ---- Derived parameters ----
omega = 2*pi*n / T;               % angular rate so that we make n turns in T
b     = (r1 - r0) / (2*pi*n);     % radial growth per radian
c     = (z1 - z0) / (2*pi*n);     % axial growth per radian

% ---- Trajectory generation ----
t     = 0:dt:T;
theta = omega * t;
rho   = r0 + b * theta;
x     = rho .* cos(theta);
y     = rho .* sin(theta);
z     = z0 + c * theta;

N=size(t,2);

f=mass*g*2;

% update limits
u_b=ones(m,1)*f/4;
% u_b=(umin+umax)/2;
umin=umin-u_b; % max q in the u_b=(umin+umax)/2, but only 0.66.
umax=umax-u_b;
plim=[umin umax];
v=[x;y;z];
% v=[x;y;z];

%% setup function of allocation lib
% ========
%% setup ACA
global NumU
NumU=m;
LPmethod=2; % LPmethod should be an integer between 0 and 5. when LPmethod=2 set upper of lambda to Inf can't save this method!!! but big number is the same as that method based linprog
% DPscaled_LPCA的结果和restoring接近但是有细微区别，对lambda限制在0-1之间，有助于优先级的理论推导。
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
%%
u=zeros(m,1);
x_LPwrap=zeros(m,N);
x_LPwrap_rest=zeros(m,N);
x_LPwrap1=zeros(m,N);
x_PCA=zeros(m,N);
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

m_higher=zeros(k,1);
Deriv_xF_v = zeros(k, k, N);
Deriv_x_v = zeros(m+1, k, N);
Deriv_u_v = zeros(m, k, N);
inF = zeros(k, N);
inD = zeros(m+1-k, N);
Deriv_u_v_inv = zeros(m, k, N);
Deriv_u_v_rest = zeros(m, k, N);
Deriv_ku_u = zeros(m, m, N);
K_u = zeros(N,1);
K_status = zeros(N,1);
elapsed_time = zeros(N,1);
%% simulate flight process  
for idx=1:N  % or x:N for debug
    
    IN_MAT(1:k,end) = v(:,idx)+m_higher; %[ 36.8125; 0;92.9776];%

    % u = LPwrap(IN_MAT); % function of ACA lib
    % u= CGIwrap(IN_MAT);
    % u = LPwrap(IN_MAT); % function of ACA lib
    % [u, ~,~,~,~,~] = DP_LPCA_copy_for_dec(v(:,idx),B,umin,umax,100);
    [u, ~,~] = DP_LPCA_prio(m_higher,v(:,idx),B,umin,umax,100);
    u=min(max(u, umin), umax);
    x_LPwrap(:,idx) =u;
    x_LPwrap_rest(:,idx) =restoring_opt(B,u,umin,umax);
    
    %% for deriv
    
    % [u, ~,Deriv_xF_v(:,:,idx),inF(:,idx),inD(:,idx),xout] = DP_LPCA_copy_for_dec(v(:,idx),B,umin,umax,100);
    % u=min(max(u, umin), umax);
    % x_LPwrap(:,idx) =u;%restoring(B,u,umin,umax);
    % tic;
    % [x_LPwrap_rest(:,idx),K_u(idx),K_status(idx),Deriv_ku_u(:,:,idx)]=restoring_opt2(B,u,umin,umax);
    % [x_LPwrap_rest(:,idx)]=restoring_opt(B,u,umin,umax);
    % elapsed_time(i) = toc;
    
    % x_LPwrap_rest(:,idx)=min(max(x_LPwrap_rest(:,idx), umin), umax);
    % 
    % Deriv_x_v(inF(:,idx)',:,idx)=Deriv_xF_v(:,:,idx);
    % Deriv_u_v(:,:,idx)=Deriv_x_v(1:m,:,idx);
    % Deriv_u_v_inv(:,:,idx)=pinv(B);
    % Deriv_u_v_rest(:,:,idx)=(eye(m,m)+ Deriv_ku_u(:,:,idx))*Deriv_u_v(:,:,idx);


    u=pinv(B)*v(:,idx);
    u=min(max(u, umin), umax);
    x_inv(:,idx) =u;% restoring(B,u,umin,umax);
    % [x_inv(:,idx),K_u(idx),K_status(idx),Deriv_ku_u(:,:,idx)] = restoring_return_k(B,u,umin,umax);

    [u,~,~] = wls_alloc(B,v(:,idx),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
    u=min(max(u, umin), umax);
    x_wls(:,idx) = u;
    % [x_wls(:,idx),K_u(idx),K_status(idx),Deriv_ku_u(:,:,idx)] = restoring_return_k(B,u,umin,umax);
    %%

    % [u, ~,~] = DP_LPCA_prio(m_higher,v(:,idx),B,umin,umax,100);
    % u=min(max(u, umin), umax);
    % x_PCA(:,idx) =restoring_cpp(B,u,umin,umax);

    % u = LPwrap(IN_MAT1); % incremental form. 
    % The control constraint δ ≤ δ ≤ δ must contain the origin, i.e., δ = 0
    % must be a feasible control input. 
    % In order word, 0 have to be a feasible solution. When add vel contraint to s.t. 
    % and if optimization variables is delta_u, then 0 is a feasible
    % solution, else if optimization variables is u, then 0 is not a
    % feasible solution, so the LP Not working.
    % x_LPwrap_incre(:,idx) = min(max(u, umin-u_0), umax-u_0)+u_0;

    % u= CGIwrap(IN_MAT);
    % u=min(max(u, umin), umax);
    % x_CGIwrapp(:,idx) = restoring(B,u,umin,umax);
    % 
    % u = DAwrap(IN_MAT);
    % u=min(max(u, umin), umax);
    % x_DAwrap(:,idx) = restoring(B,u,umin,umax);
    % 
    % u = VJAwrap(IN_MAT);
    % u=min(max(u, umin), umax);
    % x_VJAwrap(:,idx) = restoring(B,u,umin,umax);
    % 
    % u=pinv(B)*v(:,idx);
    % u=min(max(u, umin), umax);
    % x_inv(:,idx) = restoring(B,u,umin,umax);
    % 
    % [u,~,~] = wls_alloc(B,v(:,idx),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
    % u=min(max(u, umin), umax);
    % x_wls(:,idx) = restoring(B,u,umin,umax);
    % 
    % u =wls_alloc_gen(B,v(:,idx),umin,umax,eye(k),eye(m),zeros(m,1),1e6,zeros(m,1),zeros(m,1),100,4);
    % u=min(max(u, umin), umax);
    % x_wls_gen(:,idx) = restoring(B,u,umin,umax);
    % 
    % [u,~] = dir_alloc_linprog(B,v(:,idx), umin, umax, 1e4); % LPmethod=2 and lam=1 of dir_alloc_linprog is lager but similar
    % u=min(max(u, umin), umax);
    % x_dir_alloc_linprog(:,idx) = restoring(B,u,umin,umax);
    % 
    % [u,~] = dir_alloc_linprog_re(B,v(:,idx), umin, umax);
    % u=min(max(u, umin), umax);
    % x_dir_alloc_linprog_re(:,idx) = restoring(B,u,umin,umax);

    % [u,~] = dir_alloc_linprog_re_bound(B,v(:,idx), umin, umax, 1e4);% the
    % % same as dir_alloc_linprog for any lam >=1, lam have to be >1 when use
    % % linprog, that will be the same as LPmethod=3
    % u=min(max(u, umin), umax);
    % x_dir_alloc_linprog_re_bound(:,idx) = restoring(B,u,umin,umax);

    % [u,~] = use_LP_lib(B,v(:,idx), umin, umax); % ToDo: use the LP lib
    % x_use_LP_lib(:,idx)=min(max(u, umin), umax);

    % [u,~,~] =allocator_dir_LPwrap_4(single(B), single( v(:,idx)), single(umin),single(umax)); % ToDo: 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.914283e-09。 
    % x_allocator_dir_LPwrap_4(:,idx) = min(max(u, umin), umax);
    
end
fprintf('代码执行时间：%.9f 秒\n',mean( elapsed_time));


x1=x_inv; 
x2=x_LPwrap;
x3=x_LPwrap_rest;
x4=x_wls;








% actual moments produced. The B matrix have to be the same.
U1=B*x1;
U2=B*x2;
U3=B*x3;
U4=B*x4;
% figure;
% if(k==3)
% 
%     plot3(x, y, z, 'k-', 'LineWidth', 2);
%     grid on;
%     axis equal;
%     xlabel('X轴');
%     ylabel('Y轴');
%     zlabel('Z轴');
%     title('三维螺旋线');
%     view(3);    % 三维视角
% 
%     % 添加标注
%     text(x(1), y(1), z(1), '起点', 'FontSize', 12);
%     legend('螺旋线', 'Location', 'best');
% else
%     plot(x, y, 'b-', 'LineWidth', 1.5);
%     axis equal;   % 等比例坐标
%     grid on;
%     title('二维阿基米德螺旋线');
%     xlabel('X轴');
%     ylabel('Y轴');
% end
Deriv_color   = [0.10 0.35 0.75];   
inF_color   = [0.85 0.33 0.10];  
inv_color = [0.85, 0.10, 0.40]; 
d_color  = [0.00, 0.00, 0.00];      
pca_color  = [0.00, 0.60, 0.60];   
wls_color = [0.70, 0.50, 0.00];
cmd_color = [0.40, 0.00, 0.70];
 
figure,
for indx=1:m
subplot(m,1,indx)
plot(t,x1(indx,:),'--', 'Color', inv_color, 'LineWidth', 1.5);hold on;
plot(t,x2(indx,:),'-.', 'Color', d_color, 'LineWidth', 1.2);hold on;
plot(t,x3(indx,:),'-', 'Color', pca_color, 'LineWidth', 1.);hold on;
plot(t,x4(indx,:),'-','Color', wls_color,'LineWidth', 0.5);hold on;
legend('u_{inv}', 'u_{d}', 'u_{pca}','u_{wls}','Location', 'northwest','NumColumns', 4);
% legend('u_{inv}', 'u_{d}', 'u_{pca}','Location', 'northwest','NumColumns', 3);
% axis([0 2 -0.6 0.6]); xticks(0:0.4:2);yticks(-0.6:0.3:0.6);
xlabel('Time (s)');
ylabel(sprintf('u_%d', indx));
grid on;
end

figure,
for indx=1:k
subplot(k,1,indx)
plot(t,v(indx,:),':', 'Color', cmd_color, 'LineWidth', 1.8);hold on;
plot(t,U1(indx,:),'--', 'Color', inv_color, 'LineWidth', 1.4);hold on;
plot(t,U2(indx,:),'-.', 'Color', d_color, 'LineWidth', 1.8);hold on;
plot(t,U3(indx,:),'-', 'Color', pca_color, 'LineWidth', 0.5);
plot(t,U4(indx,:),'-','Color', wls_color,'LineWidth', 0.5);hold on;
legend('$\nu$','$\nu_{inv}$', '$\nu_{d}$','$\nu_{pca}$','$\nu_{wls}$','Interpreter', 'latex','Location', 'northwest','NumColumns', 1);
% legend('$\nu$','$\nu_{inv}$', '$\nu_{d}$','$\nu_{pca}$','Interpreter', 'latex','Location', 'northwest','NumColumns', 1);
% axis([0 2 -1.5 1.5]); xticks(0:0.5:2);yticks(-1.5:0.5:1.5);
xlabel('Time (s)');
if k==1
    ylabel('$\nu(1)$', 'Interpreter', 'latex');
elseif k==2
    ylabel('$\nu(2)$', 'Interpreter', 'latex');
end
grid on;
end
figure,
q=vview_norm3D(B,plim,pinv(B),v, U1,U2,U3,U4)

 