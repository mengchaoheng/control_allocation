clear all;
close all;
B=[-0.5     0       0.5     0;
                     0      -0.5     0       0.5;
                     0.25    0.25    0.25    0.25];
uMin=ones(4,1)*(-20)*pi/180;
uMax=ones(4,1)*20*pi/180;

yd=[-0.2;-0.2;0]; 
ye=[0;0;-0.2];
IN_MAT = [B yd;uMin' 0;uMax' 0];           
%===================================∑˘÷µ≤‚ ‘==================================================
N=20;
[X,Y,Z] = sphere(N);
x1=zeros(4,(N+1)^2);
u1=zeros(4,1);

for i=1:(N+1)^2%length(M_des(1:1000,1))%%length(X)
v=1*[X(i);Y(i);Z(i)];% –Èƒ‚÷∏¡ÓM_des(i,:)'%
IN_MAT(1:3,end)=v;
%=====================================
u1= LPwraparm(IN_MAT);
x1(:,i)=Constrain(u1,uMin,uMax);
end
U1=B*x1;

figure(1),
plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
