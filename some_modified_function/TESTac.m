clear all;
close all;
addpath(genpath(pwd));
global NumU Wp
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp=W'*W; %Initialize Weight matrix to unity to start
B=zeros(3,NumU);uMin=zeros(NumU,1);uMax=zeros(NumU,1);INDX=zeros(1,NumU);
effector=4;
n=3;
B(:,1:effector)=[-0.5     0       0.5     0;
                     0      -0.5     0       0.5;
                     0.25    0.25    0.25    0.25];
% B(:,1:effector)=[-1     0      1     0;
%                   0    -1      0     1;
%                   1     1      1     1];
uMin(1:effector)=ones(effector,1)*(-20)*pi/180;
uMax(1:effector)=ones(effector,1)*20*pi/180;

% B(:,1:effector) =[0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
%     1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
%    -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
% uMin(1:effector) =[-0.9599;-0.9599;-0.5236;-0.5236;-0.5236;-0.5236;-0.5236];
% uMax(1:effector) =[0.4363;0.4363;0.5236;0.5236;0.5236;0.5236;0.5236];

INDX(1:effector)=ones(1,effector);
% IN_MAT = [B     d ye;
%           umin' 0 0;
%           umax' 0 0;
%           INDX  LPmethod 0];

LPmethod=10;% just 2 3 7 ok
yd=[0;0;0]; 
ye=[0;0;0];
IN_MAT = [B yd ye;uMin' 0 0;uMax' 0 0;INDX  LPmethod 0];
[u,~,~] = LPwrap(IN_MAT); 

actuator_controls = csvread('log_23_2021-6-5-17-24-50_actuator_controls_0_0.csv',1,0);% 3 4 5
indi_feedback_input = csvread('log_23_2021-6-5-17-24-50_indi_feedback_input_0.csv',1,0);% 3 4 5
% actuator_controls = csvread('log_22_2021-6-5-17-24-16_actuator_controls_0_0.csv',1,0);% 3 4 5
% indi_feedback_input = csvread('log_22_2021-6-5-17-24-16_indi_feedback_input_0.csv',1,0);% 3 4 5
x=zeros(4,length(actuator_controls));
x1=zeros(4,length(actuator_controls));
error=zeros(3,length(actuator_controls));
error1=zeros(3,length(actuator_controls));
errord=zeros(3,length(actuator_controls));
errore=zeros(3,length(actuator_controls));
for i=1:length(actuator_controls)%%length(X)
ye=indi_feedback_input(i,3:5)';% 虚拟指令M_des(i,:)'%
yd=actuator_controls(i,3:5)'-ye;
IN_MAT(1:3,end-1)=yd;
IN_MAT(1:3,end)=ye;
%=====================================
u=pinv(B(:,1:effector))*actuator_controls(i,3:5)';
x(:,i)=Constrain(u,uMin(1:effector),uMax(1:effector));
error(:,i)=actuator_controls(i,3:5)' - B(:,1:effector)*x(:,i);


[u1,error1(:,i),errore(:,i),errord(:,i)] = LPwrap(IN_MAT);
x1(:,i)=u1(INDX>0.5);%Constrain(u1(INDX>0.5),uMin,uMax);
% [u1,~,~] = dir_alloc_sim(v, uMin(1:effector),uMax(1:effector), B(:,1:effector));
% x1(:,i)=Constrain(u1,uMin(1:effector),uMax(1:effector));
% u1=DP_LPCA(v,B(:,INDX>0.5),uMin(1:effector),uMax(1:effector),100,3,4);
% x1(:,i)=Constrain(u1,uMin(1:effector),uMax(1:effector));

% u2=SBprio_LPCA(v,ye,B,ep*ones(4,1),zeros(4,1),uMin,uMax,5e2);
% x2(:,i)=Constrain(u2,uMin,uMax);
end
U=B(:,INDX>0.5)*x;
U1=B(:,INDX>0.5)*x1;
figure,
plot(error(1,:));hold on;
plot(error(2,:));hold on;
plot(error(3,:));hold on;
figure,
plot(error1(1,:));hold on;
plot(error1(2,:));hold on;
plot(error1(3,:));hold on;
figure,
plot(errord(1,:));hold on;
plot(errord(2,:));hold on;
plot(errord(3,:));hold on;
figure,
plot(errore(1,:));hold on;
plot(errore(2,:));hold on;
plot(errore(3,:));hold on;
% figure,
% plot3(U(1,:),U(2,:),U(3,:),'b<');
% hold on;
% plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
