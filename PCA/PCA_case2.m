clear all;
close all;

addpath(genpath('/Users/mch/Proj/control_allocation'));
%% 全局图形设置
set(groot, ...
    'defaultAxesFontSize', 8, ...
    'defaultAxesFontName', 'Times New Roman', ...
    'defaultAxesLineWidth', 0.5, ...
    'defaultAxesLabelFontSizeMultiplier', 1, ...
    'defaultAxesTitleFontSizeMultiplier', 1);
%% 
B= [1    0      0.5;
     0      -1     -0.5];

[k,m] = size(B);
umin=ones(m,1)*-0.5;
umax=ones(m,1)*0.5;
plim=[umin umax];




% =========
v_h=[0.9; 0.0];
v_c=[0.0; 0.6]; 
v_c1=[-0.4; 0.4]; 
v_c2=[-1.2; 0.9]; 
% ==========
[u, ~,~] = DP_LPCA_prio(v_h,v_c,B,umin,umax,100);
u=min(max(u, umin), umax);
x_PCA =restoring_cpp(B,u,umin,umax);
v_all= B*x_PCA;

% ==========
[u, ~,~] = DP_LPCA_prio(v_h,v_c1,B,umin,umax,100);
u=min(max(u, umin), umax);
x_PCA =restoring_cpp(B,u,umin,umax);
v_all1= B*x_PCA;
% ==========
[u, ~,~] = DP_LPCA_prio(v_h,v_c2,B,umin,umax,100);
u=min(max(u, umin), umax);
x_PCA =restoring_cpp(B,u,umin,umax);
v_all2= B*x_PCA;

fig1=figure(1);

q=vview_case(B,plim,pinv(B), v_h, v_c,[], v_c1,v_all1,v_c2,v_all2)
% axis([-1 1 -1.8 1.1]);
% xticks(-1:0.5:1);
% yticks(-1.8:0.3:1.2);
axis([-1 1 -1 1]);
xticks(-1:0.5:1);
yticks(-1:0.5:1);
PlotToFileColorPDF(fig1,'results/CA_case2.pdf',6.5,6.5);




