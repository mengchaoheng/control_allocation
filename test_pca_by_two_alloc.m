clc;
clear;
close all;
% function callqpact
%===============================================================
B=[-0.5   0       0.5   0;
     0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
umin=[1;1;1;1]*(-20)*pi/180;
umax=[1;1;1;1]*20*pi/180;
%===================================幅值测试==================================================
N=40;
[XX,YY,ZZ] = sphere(N);
%==================================速度约束测试==================================
t=0:0.01:0.5;
% X1=0.3*sin(pi*t);Y1=0.3*cos(pi*t);Z1=0.1*sin(pi*t);
% X2=0.1*sin(pi*t);Y2=0.1*cos(pi*t);Z2=0.05*sin(pi*t);
X1=0.0;Y1=0.0;Z1=0.3;
X2=0.1;Y2=0.3;Z2=0;
X=X1+X2;Y=Y1+Y2;Z=Z1+Z2;
x=zeros(4,length(X1));%(N+1)^2);
xx=zeros(4,length(X1));%(N+1)^2);
x1=zeros(4,(N+1)^2);
x2=zeros(4,(N+1)^2);
x3=zeros(4,(N+1)^2);
x4=zeros(4,(N+1)^2);
iter=zeros((N+1)^2,1);
fval=zeros((N+1)^2,1);
t=zeros(3,(N+1)^2);
u=[0;0;0;0];
uu=zeros(4,1);
u1=[0;0;0;0];
u2=[0;0;0;0];
u3=[0;0;0;0];
u4=[0;0;0;0];
p_limits=20;
for i=1:(N+1)^2%
vv=1*[XX(i);YY(i);ZZ(i)];
% uu = wls_alloc_mul(vv,uu);
% xx(:,i)=uu;
% umin=max([1;1;1;1]*(-20)*pi/180,-0.01*400*pi/180+uu);
% umax=min([1;1;1;1]*20*pi/180,0.01*400*pi/180+uu);
% uu = dir_alloc_mch(vv, umin,umax);
uu = wls_alloc_mul(vv, uu,p_limits, 0);
xx(:,i)=uu;
end
for i=1:length(X1)  %(N+1)^2%
v1=1*[X1(i);Y1(i);Z1(i)]; % 虚拟指令
v2=1*[X2(i);Y2(i);Z2(i)];
v=1*[X(i);Y(i);Z(i)];
% % %==================有效集=====================
% u = dir_alloc_mch(v, umin,umax);  
u =wls_alloc_mul(v, u, p_limits, 1);
x(:,i)=u;
% u1 = dir_alloc_mch(v1, umin,umax); 
u1 =wls_alloc_mul(v1, u1, p_limits, 1);
x1(:,i)=u1;
u2 = dir_alloc_mch(v2, umin,umax);
% u2 =wls_alloc_mul(v2, u2, umin,umax);
x2(:,i)=u2;
% 改进
u3=two_dir_alloc_mch(v1, v2, p_limits,1,u3);
x3(:,i)=u3;
u4=wls_alloc_mul(v, u4, p_limits, 1);
x4(:,i)=u4;
end
UU=B*xx;
U=B*x;
U1=B*x1;
U2=B*x2;
U3=B*x3;
U4=B*x4;
% figure,
plot3(UU(1,:),UU(2,:),UU(3,:),'k.');grid on;
hold on;
plot3(X,Y,Z,'b>');
hold on;
plot3(U(1,:),U(2,:),U(3,:),'g>');grid on;
hold on;
plot3(U1(1,:),U1(2,:),U1(3,:),'r*');grid on;
% legend('边界','给定','原','改进')
hold on;
plot3(U2(1,:),U2(2,:),U2(3,:),'g+');grid on;
hold on;
plot3(U3(1,:),U3(2,:),U3(3,:),'b*');grid on;
plot3(U4(1,:),U4(2,:),U4(3,:),'r+');grid on;
hold on;
legend('边界','给定','原','扰动','动态','改进','二次')