clear all;
close all;
addpath(genpath(pwd))
folder ='some_modified_function'; 
rmpath(folder) % remove old version
folder ='s-function_used_in PlanD'; 
rmpath(folder) % remove old version
global NumU 
NumU=4; % Number of controls
INDX=ones(1,NumU);
B=[-0.5     0       0.5     0;
                     0      -0.5     0       0.5;
                     0.25    0.25    0.25    0.25];
% B(:,1:effector)=[-1     0      1     0;
%                   0    -1      0     1;
%                   1     1      1     1];
uMin=ones(4,1)*(-20)*pi/180;
uMax=ones(4,1)*20*pi/180;
use_date=0;
if(use_date)
    load 'variables.mat'; % run '/New_LP_dir/allocation_log/plot_states.m' for y_all and u_px4
    [N,~]=size(y_all);  
    x1=zeros(4,N);
    u1=zeros(4,1);
    x2=zeros(4,N);
    u2=zeros(4,1);
else
    M=50;
    x=zeros(4,(M+1)^2);
    u=zeros(4,1);
    [X,Y,Z] = sphere(M);
    x1=zeros(4,(M+1)^2);
    u1=zeros(4,1);
    x2=zeros(4,(M+1)^2);
    u2=zeros(4,1);

    N=(M+1)^2;
end
% ========
% IN_MAT = [B     d
%           umin' 0
%           umax' 0
%           INDX  0]
yd=[-0.2;-0.2;0]; 
IN_MAT = [B yd;uMin' 0;uMax' 0;INDX 0];
u1= CGIwrap(IN_MAT);
for i=1:N% (N+1)^2  for  sphere %length(M_des(1:1000,1))%%length(X)
% 
if(use_date)
    v=y_all(i,:)';
else
    v=0.34*[X(i);Y(i);Z(i)];
end
IN_MAT(1:3,end)=v;
u1= CGIwrap(IN_MAT);
x1(:,i)=Constrain(u1,uMin,uMax);
u2=pinv(B)*v;
x2(:,i)=Constrain(u2,uMin,uMax);
end
U1=B*x1;
U2=B*x2;


dt=0.01;
t=0:dt:dt*(N-1);

tt=1:1:(N-1);

if(use_date)

else
    figure,
    plot3(U1(1,:),U1(2,:),U1(3,:),'g*');
end