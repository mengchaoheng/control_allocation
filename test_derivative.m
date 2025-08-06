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

B= [1  0 -0.5;0 1  -0.5];
% l1=0.167;l2=0.069;k_v=3;
% I_x=0.01149;
% I_y=0.01153;
% I_z=0.00487;
% I=[I_x 0 0;0 I_y 0;0 0 I_z];
% %=============================4==================================
% B=I\[-l1     0       l1     0;
%      0      -l1     0       l1;
%      l2    l2    l2    l2]*k_v;
[k,m] = size(B);
umin=ones(m,1)*-0.5;
umax=ones(m,1)*0.5;
% umin=ones(m,1)*(-20)*pi/180;
% umax=ones(m,1)*20*pi/180;
plim=[umin umax];

% run Generate_input_data;
r_d=0.65;% determin by AS  10

dt=0.001;
T_len=4;
t=0:dt:T_len;
N=size(t,2);
v=zeros(k,N);
% v(1:2,:)=r_d*[sin(2*pi*t);cos(2*pi*t)];
a = 0;       % 起始半径 0
b = 0.01;     % 螺距控制参数 2
c = 0.1;    % 螺距系数 2
theta = t*10*pi;  % 角度范围（0~10π）
rho = a + b * theta;   % 极径方程
x = rho .* cos(theta); % 转直角坐标x
y = rho .* sin(theta); % 转直角坐标y
z = c * theta;  % z轴线性增长
v=[x;y];
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
tic;
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

%% simulate flight process  
for idx=1:N  % or x:N for debug
    
    IN_MAT(1:k,end) = v(:,idx)+m_higher; %[ 36.8125; 0;92.9776];%

    % u = LPwrap(IN_MAT); % function of ACA lib
    
    %% for deriv
    [u, ~,Deriv_xF_v(:,:,idx),inF(:,idx),inD(:,idx),xout] = DP_LPCA_copy_for_dec(v(:,idx),B,umin,umax,100);
    u=min(max(u, umin), umax);
    x_LPwrap(:,idx) =u;%restoring(B,u,umin,umax);
    [x_LPwrap_rest(:,idx),K_u(idx),K_status(idx),Deriv_ku_u(:,:,idx)]=restoring_return_k(B,u,umin,umax);

    

    Deriv_x_v(inF(:,idx)',:,idx)=Deriv_xF_v(:,:,idx);

    Deriv_u_v(:,:,idx)=Deriv_x_v(1:m,:,idx);
    Deriv_u_v_inv(:,:,idx)=pinv(B);
    Deriv_u_v_rest(:,:,idx)=(eye(m,m)+ Deriv_ku_u(:,:,idx))*Deriv_u_v(:,:,idx);
    if (K_u(idx)~=0 && K_status(idx)==2 )
        % Deriv_ku_u(:,:,idx)*Deriv_u_v(:,:,idx)
        % K_u(idx)
        % inF(:,idx)
        % xout(m+1)
        % disp('lambda = 1 ');
    end
    if(abs(xout(m+1)-1)<100*eps)
        % K_u(idx)
        % K_status(idx)
        % Deriv_ku_u(:,:,idx)
    end


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
elapsed_time = toc;
fprintf('代码执行时间：%.2f 秒\n', elapsed_time);


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
plot(t,x1(indx,:),'--', 'Color', inv_color, 'LineWidth', 0.8);hold on;
plot(t,x2(indx,:),'-.', 'Color', d_color, 'LineWidth', 0.8);hold on;
plot(t,x3(indx,:),'-', 'Color', pca_color, 'LineWidth', 1);hold on;
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
plot(t,v(indx,:),':', 'Color', cmd_color, 'LineWidth', 0.8);hold on;
plot(t,U1(indx,:),'--', 'Color', inv_color, 'LineWidth', 0.8);hold on;
plot(t,U2(indx,:),'-.', 'Color', d_color, 'LineWidth', 0.8);hold on;
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
q=vview_norm(B,plim,pinv(B),v, U1,U2,U3,U4)

% 画图，每个元素一个子图
%% —— 1) 数值差分重算 Deriv_u_v_rest —— 
% 假设 v 是 2×N，x3 是 3×N（即 x_LPwrap_rest），
% 并且你已经有一个函数 [u_rest] = alloc_rest(v) 能给出
% 对应的 x3(:,idx) = u_rest(v(:,idx))。

% h = 1e-6;            % 差分步长（可根据量纲微调）
% Deriv_num = zeros(m, k, N);
% for idx = 1:N
%     v0 = v(:,idx);
%     for j = 1:k
%         dv = zeros(k,1);
%         dv(j) = h;
%         % 正向、负向
%         [u, ~,~,~,~] = DP_LPCA_copy_for_dec(v0+ dv,B,umin,umax,100);
%         u=min(max(u, umin), umax);
%         [u_plus,~,~,~]=restoring_return_k(B,u,umin,umax);
%         [u, ~,~,~,~] = DP_LPCA_copy_for_dec(v0- dv,B,umin,umax,100);
%         u=min(max(u, umin), umax);
%         [u_minus,~,~,~]=restoring_return_k(B,u,umin,umax);
% 
% 
%         % 中心差分
%         Deriv_num(:,j,idx) = (u_plus - u_minus) / (2*h);
%     end
% end

%% —— 2) 绘图 —— 
% figure;
% colors_num = {[0.85 0.10 0.10], [0.10 0.30 0.80]};  % 红／蓝
% 
% for i = 1:3
%     subplot(3,1,i);
%     hold on; grid on;
% 
%     d_v1 = squeeze(Deriv_num(i,1,:));  
%     d_v2 = squeeze(Deriv_num(i,2,:));  
% 
%     plot(t, d_v1, '-',  'Color', colors_num{1}, 'LineWidth', 1.2, ...
%          'DisplayName', sprintf('\\partial x_{3,%d}/\\partial v_1', i));
%     plot(t, d_v2, '--', 'Color', colors_num{2}, 'LineWidth', 1.2, ...
%          'DisplayName', sprintf('\\partial x_{3,%d}/\\partial v_2', i));
% 
%     xlabel('Time (s)', 'Interpreter','latex');
%     ylabel('Derivative', 'Interpreter','latex');
%     title( sprintf('Component $x_{3,%d}$', i), ...
%            'Interpreter','latex', 'FontSize', 10 );
% 
%     legend('Location','best', 'Interpreter','latex','FontSize',8);
%     xlim([t(1) t(end)]);
% end
% 
% sgtitle('Numerical Jacobian $\partial x_3/\partial v$', ...
%         'Interpreter','latex','FontSize',12);

figure;
for i = 1:m
    for j = 1:k
        subplot(m*k, 1, (i-1)*k + j);
        plot(t, squeeze(Deriv_u_v(i,j,:)), '-','Color', Deriv_color, 'LineWidth', 0.8);hold on;
        plot(t, squeeze(Deriv_u_v_inv(i,j,:)), '--','Color', inv_color, 'LineWidth', 0.8);hold on;
        plot(t, squeeze(Deriv_u_v_rest(i,j,:)), '-.','Color', d_color, 'LineWidth', 1.2);hold on;
        % plot(t, squeeze(Deriv_num(i,j,:)), 'g:', 'LineWidth', 2);hold on;
        plot(t, K_u*10, 'm-', 'LineWidth', 0.5);hold on;
        plot(t, K_status, 'c-', 'LineWidth', 0.5);hold on;
        xlabel('Time (s)');
        ylabel(sprintf('du/dv(%d,%d)', i, j));
        legend('Deriv_{LP}', 'Deriv_{inv}', 'Deriv_{rest}','k_u*10', 'k_status');
        axis([1 3 -inf inf]);
        grid on;
        
    end
end
% 画图，每个元素一个子图
figure;
for i = 1:k
        subplot(k, 1, i);
        plot(t, squeeze(inF(i,:)), '-','Color', inF_color, 'LineWidth', 0.8);hold on;
        xlabel('Time (s)');
        ylabel(sprintf('index %d of F', i));
        grid on;
end

figure;
plot(t, K_u, '-','Color', inF_color, 'LineWidth', 0.8);hold on;
plot(t, K_status, '-','Color', d_color, 'LineWidth', 0.8);hold on;
axis([1.8 3 -inf inf]);