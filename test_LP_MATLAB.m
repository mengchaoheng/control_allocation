clear all;
close all;
addpath(genpath(pwd))
folder ='some_modified_function'; 
rmpath(folder) % remove old version
folder ='s-function_used_in PlanD'; 
rmpath(folder) % remove old version
% 注意B不同
% B=[-0.5     0       0.5     0;
%      0      -0.5     0       0.5;
%      0.25    0.25    0.25    0.25];

l1=0.148;l2=0.069;k_v=3;
B=k_v*[-l1     0       l1     0;
     0      -l1     0       l1;
     l2    l2    l2    l2];

B_inv=pinv(B)
[k,m] = size(B);
% m=4;
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
use_date=0; % test use px4 fly data.
if(use_date)
    % load 'hover.mat'; % run '/New_LP_dir/allocation_log/plot_states.m' for y_all and u_px4
    load 'fly.mat';
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


for i=1:N% (N+1)^2  for  sphere %length(M_des(1:1000,1))%%length(X)
% 
    if(use_date)
        v=(y_all(i,:)');
    else
        v=(0.2*[X(i);Y(i);Z(i)]);
    end
    [u1,~] = dir_alloc_linprog_re(B,v, umin, umax);
    x1(:,i) = Constrain(u1,umin,umax);


    [u2,~] = test_LP_lib(B,v, umin, umax);
    x2(:,i)=Constrain(u2,umin,umax);
end
U1=B*x1;
U2=B*x2;



dt=0.01;
t=0:dt:dt*(N-1);

tt=1:1:(N-1);
if(use_date)
    y_all_new=y_all';
    U_px4=B*u_px4';
    error1=U1-y_all_new;
    error2=U2-y_all_new;
    error3=U_px4-y_all_new;
    figure,
    subplot(4,1,1)
    plot(t,x1(1,:),'r.');hold on;
    plot(t,x2(1,:),'b--');hold on;
    % plot(t,u_px4(:,1),'g-');hold on;
    subplot(4,1,2)
    plot(t,x1(2,:),'r.');hold on;
    plot(t,x2(2,:),'b--');hold on;
    % plot(t,u_px4(:,2),'g-');hold on;
    subplot(4,1,3)
    plot(t,x1(3,:),'r.');hold on;
    plot(t,x2(3,:),'b--');hold on;
    % plot(t,u_px4(:,3),'g-');hold on;
    subplot(4,1,4)
    plot(t,x1(4,:),'r.');hold on;
    plot(t,x2(4,:),'b--');hold on;
    % plot(t,u_px4(:,4),'g-');hold on;
    figure,
    subplot(3,1,1)
    plot(t,error1(1,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,error2(1,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,error1(1,:)-error2(1,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    % plot(t,U1(1,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,U2(1,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,y_all(:,1),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;

    subplot(3,1,2)
    plot(t,error1(2,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,error2(2,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,error1(2,:)-error2(2,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    % plot(t,U1(2,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,U2(2,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,y_all(:,2),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    subplot(3,1,3)
    plot(t,error1(3,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,error2(3,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,error1(3,:)-error2(3,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    % plot(t,U1(3,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,U2(3,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,y_all(:,3),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
else
    figure,
    plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
    figure,
    plot3(U2(1,:),U2(2,:),U2(3,:),'g*');
end