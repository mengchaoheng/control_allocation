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
% setup LPwrap
% global NumU 
NumU=m; % Number of controls
LPmethod=3; % LPmethod should be an integer between 0 and 5
INDX=ones(1,m);  % active effectors
IN_MAT = [B     zeros(k,1)
          umin' 0
          umax' 0
          INDX  LPmethod];
% test for LPwrap
u = LPwrap(IN_MAT);      
u(INDX>0.5)

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

% [u1,~] = dir_alloc(B,v,umin,umax); % function of qcat lib
IN_MAT(1:3,end)=v; u1 = LPwrap(IN_MAT); % function of ACA lib

% u1 = dir_alloc_simplex(B, v, umin,umax, m); % -- mch
% u1=dir_ca_4df(v, umin,umax);  % for ductedfan whit 4 control surface -- mch 
% [u1,~] = dir_linprog_ca_4df(B,v, umin, umax); % for 4df by matlat linprog function --mch

x1(:,i)=Constrain(u1,umin,umax);

end
U1=B*x1;
% U=B*x;
plot3(U1(1,:),U1(2,:),U1(3,:),'g*');

