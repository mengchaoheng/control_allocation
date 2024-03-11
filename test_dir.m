clear all;
close all;
folder ='some_modified_function';  % 
addpath( genpath(folder) );

B=zeros(3,4);uMin=zeros(4,1);uMax=zeros(4,1);INDX=zeros(1,4);

B=[-0.5     0       0.5     0;
                     0      -0.5     0       0.5;
                     0.25    0.25    0.25    0.25];
% B(:,1:effector)=[-1     0      1     0;
%                   0    -1      0     1;
%                   1     1      1     1];
uMin=ones(4,1)*(-20)*pi/180;
uMax=ones(4,1)*20*pi/180;


      
%===================================��ֵ����==================================================
N=50;
x=zeros(4,(N+1)^2);
u=zeros(4,1);
[X,Y,Z] = sphere(N);

x1=zeros(4,(N+1)^2);
u1=zeros(4,1);


for i=1:(N+1)^2%length(M_des(1:1000,1))%%length(X)
v=5*[X(i);Y(i);Z(i)];% ����ָ��M_des(i,:)'%

%=====================================
% u=pinv(B)*v;
% x(:,i)=Constrain(u,uMin,uMax);

% [u1] = LPwrap(IN_MAT);
% x1(:,i)=u1(INDX>0.5);%Constrain(u1(INDX>0.5),uMin,uMax);
[u1,~,~] = dir_alloc_sim(v, uMin,uMax, B);
x1(:,i)=Constrain(u1,uMin,uMax);
% u1=DP_LPCA(v,B(:,INDX>0.5),uMin(1:effector),uMax(1:effector),100,3,4);
% x1(:,i)=Constrain(u1,uMin(1:effector),uMax(1:effector));

% u2=SBprio_LPCA(v,ye,B,ep*ones(4,1),zeros(4,1),uMin,uMax,5e2);
% x2(:,i)=Constrain(u2,uMin,uMax);
end
% U=B(:,INDX>0.5)*x;
U1=B*x1;



% U2=B*x2;
% V=(yd+ye);
% figure(1),
% plot3(U(1,:),U(2,:),U(3,:),'b*');
% hold on;
plot3(U1(1,:),U1(2,:),U1(3,:),'g*');
% hold on;
% plot3(ye(1,1),ye(2,1),ye(3,1),'b*');
% hold on;
% plot3(V(1,1),V(2,1),V(3,1),'g*');
% hold on;
% plot3(vv(1,:),vv(2,:),vv(3,:),'k>');

% hold on;
% plot3(u2(1,:),u2(2,:),u2(3,:),'effector*');
% 
% hold on;
% plot3(u4(1,:),u4(2,:),u4(3,:),'g>');
