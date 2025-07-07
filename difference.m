clear all;
close all;
addpath(genpath(pwd))

folder ='QCAT/qcat';  % vview所在文件夹
addpath( genpath(folder) ); 

 
B= [1    0      -0.5;
     0      1     -0.5];

[k,m] = size(B);
umin=ones(m,1)*-0.4;
umax=ones(m,1)*0.4;
plim=[umin umax];

v_h=[0.5; 0];
v_c=[0.4; 0.4];
v=v_h+v_c;


% =========
u=pinv(B)*v;
u=min(max(u, umin), umax);
x_inv = restoring(B,u,umin,umax);
v_inv= B*x_inv;
% ===========
u =wls_alloc_gen(B,v,umin,umax,eye(k),eye(m),zeros(m,1),1e6,zeros(m,1),zeros(m,1),100,3);
u=min(max(u, umin), umax);
x_wls_gen = restoring(B,u,umin,umax);
v_wls=B*x_wls_gen;
% ==========
[u, ~,~] = DP_LPCA_prio([0; 0],v,B,umin,umax,100);
u=min(max(u, umin), umax);
x_dir =restoring_cpp(B,u,umin,umax);
v_dir= B*x_dir;
% ==========
[u, ~,~] = DP_LPCA_prio(v_h,v_c,B,umin,umax,100);
u=min(max(u, umin), umax);
x_PCA =restoring_cpp(B,u,umin,umax);
v_pca= B*x_PCA;



q=vview2(B,plim,pinv(B), v_h, v_c, v, v_inv,v_wls,v_dir,v_pca)




