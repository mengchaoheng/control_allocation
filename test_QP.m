clear all;
close all;
addpath(genpath(pwd))
folder ='some_modified_function'; 
rmpath(folder) % remove old version
folder ='s-function_used_in PlanD'; 
rmpath(folder) % remove old version

B=[-0.5   0       0.5   0;
    0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
[k,m] = size(B);
% m=4;
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
%===================================幅值测试==================================================
N=50;
[X,Y,Z] = sphere(N);
%==================================速度约束测试==================================
t=0:0.01:1;
% X=0.33*sin(pi*t);
% Y=0.34*cos(pi*t);
% Z=0.17*sin(pi*t);
x=zeros(m,(N+1)^2);
x1=zeros(m,(N+1)^2);
iter=zeros((N+1)^2,1);
fval=zeros((N+1)^2,1);
t=zeros(3,(N+1)^2);
u=zeros(m,1);
u1=zeros(m,1);
for i=1:(N+1)^2%length(X)
v=0.5*[X(i);Y(i);Z(i)];% 虚拟指令
%==================有效集=====================
% test core function
% [u,~,~] =wls_alloc(B,v,umin,umax,eye(3),eye(4),zeros(4,1),1e6,zeros(4,1),zeros(4,1),100);
% u =wls_alloc_gen(B,v,umin,umax,eye(3),eye(4),zeros(4,1),1e6,zeros(4,1),zeros(4,1),100,4);
% test old version, maybe use for a simulink block
% u= qp_ca_mch([v;u],B,[umin umax],[],0.01,eye(3),eye(4),zeros(4,1),100,1e6,true);
% test 4df allocation
 u = wls_ca_4df_pv_limit(v, u, 20, 0);

% u = wls_ca_4df(v, u);
% u=dyn_ca_4df(v,u);


x(:,i)=u;
u1=pinv(B)*v;
x1(:,i)=Constrain(u1,umin,umax);
end
U=B*x;
U1=B*x1;
figure(4),
plot3(U(1,:),U(2,:),U(3,:),'b*');
hold on;
% plot3(U1(1,:),U1(2,:),U1(3,:),'r*');




% plot3(X,Y,Z,'b>');%zeros(length(X),1)
% hold on;
% plot3(UU(1,:),UU(2,:),UU(3,:),'b*');
% C=sqrt(U1.*U1+U2.*U2+U3.*U3);
% figure,
% surf(U1,U2,U3,C);% 'FaceColor','b','FaceAlpha',0.5,,'EdgeColor','none'