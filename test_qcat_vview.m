clear all;
close all;
addpath(genpath(pwd))


% 
% 1）添加当前文件夹的路径
% addpath(pwd)
% 2）删除当前文件夹的路径
% rmpath(pwd)
% 3）添加当前文件夹以及所有子文件夹的路径
% addpath(genpath(pwd))
% 4）移除当前文件夹以及所有子文件夹的路径
% rmpath(genpath(pwd))

folder ='QCAT/qcat';  % vview所在文件夹
addpath( genpath(folder) ); 

% l1=0.149;l2=0.0698;k_v=3; % k_v*delta=F on cs
% I_x=0.00967;I_y=0.0097;I_z=0.00448;
% I=diag([I_x;I_y;I_z]);
%=============================6==================================
% d=60*pi/180;
% B=I\[-l1 -l1*cos(d) l1*cos(d) l1 l1*cos(d) -l1*cos(d);0 -l1*sin(d) -l1*sin(d) 0 l1*sin(d) l1*sin(d);l2 l2 l2 l2 l2 l2]*k_v;
% B_inv=[-1 0 1;-1 -1 1;1 -1 1;1 0 1;1 1 1;-1 1 1];B=pinv(B_inv);
% 
%=============================4==================================
% B=I\[-l1 0 l1 0;0 -l1 0 l1;l2 l2 l2 l2]*k_v; 
% B=I\diag([2*l1;2*l1;4*l2])*k_v*[-0.5 0 0.5 0;0 -0.5 0 0.5;0.25 0.25 0.25 0.25];

B=[-0.5 0 0.5 0;0 -0.5 0 0.5;0.25 0.25 0.25 0.25];%= piv([-1 0 1;0 -1 1;1 0 1;0 1 1])
% I\[2*l1 0 0;0 2*l1 0;0 0 4*l2]*k_v is the different of gain, that is diag([92.4509;92.1649;186.9643])
% B=[-l1*k_v/I_x 0 l1*k_v/I_x 0;0 -l1*k_v/I_y 0 l1*k_v/I_y;l2*k_v/I_z l2*k_v/I_z l2*k_v/I_z l2*k_v/I_z];
%============================8===================================


%================================================================
[k,m] = size(B);
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
plim=[umin umax];

q=vview(B,plim,pinv(B))
% q=vview(B,plim,B_inv)


% l1=0.148;l2=0.069;k_v=3;01
% B=k_v*[-l1     0       l1     0;
%      0      -l1     0       l1;
%      l2    l2    l2    l2];
% [k,m] = size(B);
% 
% umin=ones(m,1)*(-20)*pi/180;
% umax=ones(m,1)*20*pi/180;
% plim=[umin umax];
% q1=vview(B,plim,pinv(B))

