clear all;
close all;
addpath(genpath(pwd))
folder ='some_modified_function'; 
rmpath(folder) % remove old version
folder ='s-function_used_in PlanD'; 
rmpath(folder) % remove old version



B=[-0.5     0       0.5     0;
                     0      -0.5     0       0.5;
                     0.25    0.25    0.25    0.25];
% B(:,1:effector)=[-1     0      1     0;
%                   0    -1      0     1;
%                   1     1      1     1];

[k,m] = size(B);
% m=4;
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
      
N=50;
x=zeros(m,(N+1)^2);
u=zeros(m,1);
[X,Y,Z] = sphere(N);

x1=zeros(m,(N+1)^2);
u1=zeros(m,1);


for i=1:(N+1)^2%length(M_des(1:1000,1))%%length(X)
v=0.5*[X(i);Y(i);Z(i)];% 

%=====================================
u=pinv(B)*v;
x(:,i)=Constrain(u,umin,umax);

u1 = dir_alloc_simplex(B, v, umin,umax, m);
x1(:,i)=Constrain(u1,umin,umax);

end
U1=B*x1;
% U=B*x;
plot3(U1(1,:),U1(2,:),U1(3,:),'g*');

