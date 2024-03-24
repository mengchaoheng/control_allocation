clc;
clear all;
close all;
% If we use the Standard Forms for Linear Programming Problems
% min c'x subj. to A*x <= b
%                    0 <= x
%% so we have to reformula the direction-preserving control allocation
% problem to:
% min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
%                                umin <= u <= umax
%                                   0 <= a
% and set x=u-umin, then
% min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
%                                        x <= umax-umin
%                                        0 <= x 
%                                        0 <= a
% min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] <= -B*umin
%                              [-B v][x; a] <=  B*umin
%                              [I  0][x; a] <= umax-umin
%                                        -x <= 0 
%                                        -a <= 0
% set X=[x;a]
l1=0.148;l2=0.069;k_v=3;
B=k_v*[-l1     0       l1     0;
     0      -l1     0       l1;
     l2    l2    l2    l2];
[k,m] = size(B);
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
v=[0.2;0.1;0.1];
A=[B -v;-B v; eye(m) zeros(m,1)];
b=[-B*umin; B*umin; umax-umin];
c=[zeros(m,1); -1];
[X,fval,exitflag,output,lambda]= linprog(c,A,b,[],[],zeros(1+m,1),[]);
u = X(1:m)+umin;
a = X(m+1);
% Scale down u if a>1
if a>1
    u = u/a;
end
u

%% for SDLP 
% A b c is 
% A=[B -v;-B v; eye(m) zeros(m,1);-eye(m) zeros(m,1);zeros(1,m) -1]
% b=[-B*umin; B*umin; umax-umin;0;0;0;0;0]
% c=[zeros(m,1); -1]
% or
A
b
c

% setup A b and c to sdlp_example and then run :
% ➜  build git:(main) ✗ make          
% [ 50%] Building CXX object CMakeFiles/sdlp_example.dir/example/sdlp_example.cpp.o
% [100%] Linking CXX executable sdlp_example
% [100%] Built target sdlp_example
% ➜  build git:(main) ✗ ./sdlp_example
% optimal sol: 0.155327 0.426713   0.6981   0.6981  1.20496
% u: -0.193773 0.0776133     0.349     0.349
% optimal obj: -1.20496
u_SDLP=[-0.193773; 0.0776133;     0.349;     0.349];
% optimal obj is the same, but B*u_SDLP != v, u != u_SDLP   
B*u_SDLP


% for tiny_LP
% we use the Standard Forms for Linear Programming Problems
% min c'x subj. to A*x =b
%                  0 <= x
%% so we have to reformula the direction-preserving control allocation
% problem to:
% min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
%                                umin <= u <= umax
%                                   0 <= a
% and set x=u-umin, then
% min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
%                                0 <= x <= umax-umin
%                                0 <= a
% add slack to converted inequalities into equalities of x:
% min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
%                                     x + s = umax-umin
%                                         0 <= x 
%                                         0 <= a
%                                         0 <= s
% set X=[x; a; s], that is:
% min z=[0; -1; 0]'[x; a; s]   s.t.  [B -v 0; I 0 I][x; a; s] = [-B*umin; umax-umin] 
%                                         0 <= x 
%                                         0 <= a
%                                         0 <= s
% A=[B -v 0; I 0 I]; b=[-B*umin; umax-umin]; c=[0; -1; 0];
[k,m] = size(B);
Aeq=[B -v zeros(k,m); eye(m) zeros(m,1) eye(m)]; 
beq=[-B*umin; umax-umin]; 
c=[zeros(m,1); -1; zeros(m,1)];
% build and run the test for lp_tiny lib
% $ gcc -o test_allocation.out lp_tiny.c test_allocation.c -llapack -lblas  
% $ ./test_allocation.out
% and run more...
use_date=1; % test use px4 fly data.
if(use_date)
    % load 'hover.mat'; % run '/New_LP_dir/allocation_log/plot_states.m' for y_all and u_px4
    load 'fly.mat';
    [N,~]=size(y_all);  
    x1=zeros(4,N);
    u1=zeros(4,1);
    x2=zeros(4,N);
    u2=zeros(4,1);
    v=zeros(3,N);

else
    M=50;
    x=zeros(4,(M+1)^2);
    u=zeros(4,1);
    [X,Y,Z] = sphere(M);
    x1=zeros(4,(M+1)^2);
    u1=zeros(4,1);
    x2=zeros(4,(M+1)^2);
    u2=zeros(4,1);
    v=zeros(3,(M+1)^2);

    N=(M+1)^2;
end


for i=1:N% (N+1)^2  for  sphere %length(M_des(1:1000,1))%%length(X)
% 
    if(use_date)
        v(:,i)=(y_all(i,:)');
    else
        v(:,i)=(0.5*[X(i);Y(i);Z(i)]);
    end
    [u1,~] = dir_alloc_linprog_re(B,v(:,i), umin, umax);
    x1(:,i) = Constrain(u1,umin,umax);


    % [u2,~] = test_LP_lib(B,v(:,i), umin, umax);
    % x2(:,i)=Constrain(u2,umin,umax);
end
U1=B*x1;

filename = 'data.csv';
writematrix(v',filename);
% copy the data.csv to root folder of lp_tiny and run test_input_alloc
% get output.csv to here
output = readmatrix('output.csv');
x2=output';
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