clc;
clear;
% ��������Ϊ�²�����
p_limits=[20*pi/180 4.2513];%[rad  rad/s]% 2*k_TS*speed*200=4.2513 % 200Ϊ�����ṩ��ת���ص�ת�ٲ�+-200
v_limits=[400*pi/180 9.9197];% [rad/s rad/s^2]% ʵ��ת�����仯��466.6667rad/s^2��2*k_TS*speed*466.6667=9.9197
%===============================================================
k_TS=9.9796018325697625989171178675552e-6;
speed=1065;
kc=3.157;
l_1=0.17078793-0.09;% roll,pitch
l_2=0.175;% T
l_3=0.06647954;% yaw1 yaw3
l_4=0.06647954+0.175;% yaw4
l_5=0.175-0.06647954;% yaw2
I_x=0.054593;
I_y=0.017045;
I_z=0.049226;
I=[I_x 0 0;0 I_y 0;0 0 I_z];
% =========================10 sureface=================================
L=[-l_1 0 l_1 0 -l_1 0 l_1 0 l_2;
    0 -l_1 0 l_1 0 -l_1 0 l_1 0;
    l_3 -l_5 l_3 l_4 l_3 l_4 l_3 -l_5 0];
F=diag([kc kc kc kc kc kc kc kc 1],0);%
B=I\L*F;
B=[-4.6718         0    4.6718         0   -4.6718         0    4.6718         0    3.2055;
         0  -14.9632         0   14.9632         0  -14.9632         0   14.9632         0;
    4.2635   -6.9597    4.2635   15.4868    4.2635   15.4868    4.2635   -6.9597         0];
%========================================
[m,k] = size(B);
% % % % % W2_c=2*p_limits(1)/v_limits(1);
% % % % % W2_t=2*p_limits(2)/v_limits(2);
% % % % % dca_W2=1*diag([W2_c W2_c W2_c W2_c W2_c W2_c W2_c W2_c W2_t]);
% % % % % dca_W1=1*eye(9);% 
% % % % % S=zeros(k,m);
% % % % % P_1=pinv([B(:,2) B(:,4) B(:,6) B(:,8) B(:,9)]);
% % % % % S([2 4 6 8 9],:)=P_1;
% % % % %========================================
umin=[[1;1;1;1;1;1;1;1]*(-p_limits(1));-p_limits(2)];
umax=[[1;1;1;1;1;1;1;1]*p_limits(1);p_limits(2)];
vmin=[[1;1;1;1;1;1;1;1]*(-v_limits(1));-v_limits(2)];
vmax=[[1;1;1;1;1;1;1;1]*v_limits(1);v_limits(2)];
plim=[umin umax];
rlim=[vmin vmax];
gam=1e6;
Wv=eye(3);
Wu=eye(m); 
imax=100;
W2_c=2*p_limits(1)/v_limits(1);
W2_t=2*p_limits(2)/v_limits(2);
% dca_W1=diag([W1_c W1_c W1_c W1_c W1_c W1_c W1_c W1_c W1_t]);% eye(9);%
dca_W2=1*diag([W2_c W2_c W2_c W2_c W2_c W2_c W2_c W2_c W2_t]);
dca_W1=1*eye(9);% 
S=zeros(k,m);
P_1=pinv([B(:,2) B(:,4) B(:,6) B(:,8) B(:,9)])
% P_1=pinv([B(:,4) B(:,6) B(:,9)]);
% P_1=pinv(B(:,1:8));
% S=pinv(B)
S([2 4 6 8 9],:)=P_1;
% p_limits=20*pi/180;%[rad  rad/s]
% B=[-4.6718         0    4.6718         0   -4.6718         0    4.6718         0    ;
%          0  -14.9632         0   14.9632         0  -14.9632         0   14.9632         ;
%     4.2635   -6.9597    4.2635   15.4868    4.2635   15.4868    4.2635   -6.9597         ];
% %========================================
% %========================================
% umin=[1;1;1;1;1;1;1;1]*(-p_limits);
% umax=[1;1;1;1;1;1;1;1]*p_limits;

% B=[-0.5   0       0.5   0;
%     0  -0.5    0       0.5;
%     0.25   0.25   0.25   0.25];
% umin=[1;1;1;1]*(-20)*pi/180;
% umax=[1;1;1;1]*20*pi/180;
% [m,k] = size(B);
%========================================
%===================================��ֵ����==================================================
N=100;
[X,Y,Z] = sphere(N);
%==================================�ٶ�Լ������==================================
t=0:0.01:1;
% X=0.33*sin(pi*t);
% Y=0.34*cos(pi*t);
% Z=0.17*sin(pi*t);
x=zeros(k,(N+1)^2);
x1=zeros(k,(N+1)^2);
iter=zeros((N+1)^2,1);
fval=zeros((N+1)^2,1);
t=zeros(3,(N+1)^2);
u=zeros(k,1);
u1=zeros(k,1);
for i=1:(N+1)^2%length(X)
v=50*[X(i);Y(i);Z(i)];% ����ָ��
%==================��Ч��=====================
% [u] = wls_alloc_mch(v, u, p_limits, false);
%  [u] =wls_alloc_m(v, u,false);
%  [u,~,~] = wls_alloc(B,v,umin,umax);
% [u,~,~] =wls_alloc(B,v,umin,umax,eye(3),eye(9),zeros(9,1),1e6,zeros(9,1),zeros(9,1),100);
% u= qp_ca_mch([v;u],B,plim,rlim,0.01,eye(3),eye(9),zeros(9,1),100,1e6,true);
[u] = wls_alloc_m(B,v,u,p_limits,v_limits,0.01,eye(3),eye(9),zeros(9,1),100,1e6,true);
% u=dyn_alloc_m(B,v,u,p_limits,v_limits,0.01,eye(3),dca_W1,dca_W2,S,imax,gam,true);
% u = dyn_alloc_m(v,u,false);
x(:,i)=u;
% [u1] = wls_alloc_mch(v, u1, p_limits, true);
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