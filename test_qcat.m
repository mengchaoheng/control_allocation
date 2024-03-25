clc;
clear all;
close all;
addpath(genpath(pwd))


% setup 4df
B=[-0.5     0       0.5     0;
     0      -0.5     0       0.5;
     0.25    0.25    0.25    0.25];

[k,m] = size(B);
% m=4;
plim=[];% clear before
plim(:,1)=ones(m,1)*(-20)*pi/180;
plim(:,2)=ones(m,1)*20*pi/180;
rlim=[];
% setup input data
use_hover_date=0;
if(use_hover_date)
    load 'hover.mat'; % run '/New_LP_dir/allocation_log/plot_states.m' for y_all and u_px4
    [N,~]=size(y_all);  
    u=zeros(4,N);
    u1=zeros(4,N);
    A = zeros(m,N);
    iter = zeros(1,N);
    time = zeros(1,N);
    v=zeros(k,N);
    v=y_all(:,:)';
else
    M=20;
    [X,Y,Z] = sphere(M);
    u=zeros(4,(M+1)^2);
    u1=zeros(4,(M+1)^2);
     

    N=(M+1)^2;
    v=zeros(k,N);
    for i=1:N
        v(:,i)=0.4*[X(i);Y(i);Z(i)];
    end
    A = zeros(m,N);
    iter = zeros(1,N);
    time = zeros(1,N);
end

% setup alloc_sim, see alloc_sim.m
% method   -control allocation method: 'qp'   l2-optimal allocation
%				                       'dyn'  dynamic allocation
%				                       'dir'  direct allocation
method='qp';
% Set default values of optional arguments
switch method
    case {'qp','dyn'}
% 'alg'    -numerical algorithm: 'sls'    SLS_ALLOC
%                                'wls'    WLS_ALLOC (default)
%                                'wlsc'   WLSC_ALLOC (need lib)
%                                'mls'    MLS_ALLOC
%                                'ip'     IP_ALLOC
%                                'cgi'    CGI_ALLOC
%                                'fxp'    FXP_ALLOC
        alg  = 'wls';
    case 'dir'
        alg = 'dir';
end

% Coplanar controls? Only in SLS_ALLOC
copl = iscoplanar(B);
imax = 100;	     % no of iterations
gam  = 1e6;	     % weight
tol  = 1e-6;       % tolerance in IP solver
hs   = 1;	     % hotstart
ui   = [];	     % initial control
Wi   = zeros(m,1); % initial working set
Wv   = eye(k);     % QP allocation
Wu   = eye(m);
ud   = zeros(m,1);
W1   = eye(m);     % Dynamic allocation
W2   = zeros(m);
S    = pinv(B);

if strcmp(alg,'ip')
    % Weighting matrices must be unit matrices.
    if any(any(Wv ~= eye(k))) | any(any(Wu ~= eye(m)))
      disp(' ')
      disp(['** Warning: Non-unit matrices Wv and Wu not handled by' ...
        ' IP_ALLOC **']) 
      disp(' ')
    end
    % Dynamic allocation almost always requires non-unit matrix Wu.
    if strcmp(method,'dyn')
      disp(' ')
      disp(['** Warning: Dynamic allocation not handled by IP_ALLOC,' ...
        ' resorting to WLS_ALLOC **'])
      disp(' ');
      alg = 'wls';
    end
end

if isempty(ui)
    % Set default initial conditions
    switch method
         case 'qp'
          [ui,Wi] = wls_alloc(B,v(:,1),plim(:,1),plim(:,2),Wv,Wu,ud);
         case 'dyn'
          [ui,Wi] = wls_alloc(B,v(:,1),plim(:,1),plim(:,2),Wv,W1,S*v(:,1));
         case 'dir'
          ui = dir_alloc(B,v(:,1),plim(:,1),plim(:,2));
    end
end

% Hotstart?
if hs
    u0 = ui;
    W0 = Wi;
end

% Precompute matrices used in dynamic allocation
if strcmp(method,'dyn')
    W1sq = W1^2;
    W2sq = W2^2;
    Wu = sqrtm(W1sq+W2sq);
    invWusq = inv(W1sq+W2sq);
end

% Allocate control signals to produce the demanded virtual
% control trajectory.
uprev = ui;
 
% Allocate control signals to produce the demanded virtual control trajectory.
% ========
for i=1:N% (N+1)^2  for  sphere %length(M_des(1:1000,1))%%length(X)
    
    % Compute feasible upper and lower bounds.
    if isempty(rlim)
      umin = plim(:,1);
      umax = plim(:,2);
    else
      umin = max(plim(:,1),uprev+rlim(:,1)*T);
      umax = min(plim(:,2),uprev+rlim(:,2)*T);
    end
    if ~hs
      u0 = (umin+umax)/2;
      W0 = zeros(m,1);
    end
    % Update u0 to reflect working set. Crucial when rate limits
    % are included and hotstart is used.
    i_min = W0 == -1;
    i_max = W0 == +1;
    u0(i_min) = umin(i_min);
    u0(i_max) = umax(i_max);
    % For dynamic allocation, merge the position and rate terms
    % into one term ||Wu(u-ud)||.
    if strcmp(method,'dyn')
      us = S*v(:,i);
      ud = invWusq*(W1sq*us+W2sq*uprev);
    end
    tic;
    %============================qcat===================================
    switch alg
        case 'sls'
            [u(:,i),W,iter(i)] = sls_alloc(B+copl*j,v(:,i),umin,umax,Wv,Wu,ud,u0,W0,imax);
        case 'mls'
            [u(:,i),W,iter(i)] = mls_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,u0,W0,imax);
        case 'wls'
            [u(:,i),W,iter(i)] = wls_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
        case 'wlsc'
            [u(:,i),W,iter(i)] = wlsc_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,W0,imax); % can't test by now
        case 'ip'
            [u(:,i),iter(i)] = ip_alloc(B,v(:,i),umin,umax,ud,gam,tol,imax);
        case 'fxp'
            u(:,i) = fxp_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,gam,u0,imax);
        case 'cgi'
            u(:,i) = cgi_alloc(B,v(:,i),umin,umax,Wv,Wu,ud,imax);
        case 'dir'
            u(:,i) = dir_alloc(B,v(:,i),plim(:,1),plim(:,2));
        otherwise
            error(sprintf('Unknown allocation algorithm: %s',alg));
    end
    % Register elapsed time
    time(i) = toc;
    % Limit control (only necessary with rate limits and direct alloc)
    u(:,i) = max(umin,min(umax,u(:,i)));
    % Determine active constraints in the final point.
    % +/- : max or min
    % 1/2 : position or rate limit
    A(:,i) = - 2*(u(:,i)==umin) + (u(:,i)==plim(:,1)) ...
	     + 2*(u(:,i)==umax) - (u(:,i)==plim(:,2));
    % Update uprev
    uprev = u(:,i);
    % Hotstart?
    if hs
      u0 = u(:,i);
      if any(strcmp(alg,{'sls','wls','wls_c','wlsc','mls'}))
	      W0 = W;
      end
    end
    %============================qcat===================================
    u1(:,i)  = pinv(B)*v(:,i);
    u1(:,i) = max(umin,min(umax,u1(:,i)));
end
%============================inv===================================
% moment by effector
U=B*u;
U1=B*u1;
% plot
dt=0.01;
t=0:dt:dt*(N-1);
tt=1:1:(N-1);
if(use_hover_date)
    error=U-v; % allocation error
    error1=U1-v;
    figure,
    subplot(4,1,1)
    plot(t,u(1,:),'r-');hold on;
    plot(t,u1(1,:),'b--');hold on;
    % plot(t,u_px4(:,1),'g.');hold on;
    subplot(4,1,2)
    plot(t,u(2,:),'r-');hold on;
    plot(t,u1(2,:),'b--');hold on;
    % plot(t,u_px4(:,2),'g.');hold on;
    subplot(4,1,3)
    plot(t,u(3,:),'r-');hold on;
    plot(t,u1(3,:),'b--');hold on;
    % plot(t,u_px4(:,3),'g.');hold on;
    subplot(4,1,4)
    plot(t,u(4,:),'r-');hold on;
    plot(t,u1(4,:),'b--');hold on;
    % plot(t,u_px4(:,4),'g.');hold on;
    figure,
    subplot(3,1,1)
    plot(t,error(1,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
    plot(t,error1(1,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,error(1,:)-error1(1,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    % plot(t,U1(1,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,U2(1,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,y_all(:,1),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;

    subplot(3,1,2)
    plot(t,error(2,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
    plot(t,error1(2,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,error(2,:)-error1(2,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    % plot(t,U1(2,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,U2(2,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,y_all(:,2),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    subplot(3,1,3)
    plot(t,error(3,:),'Color','r','LineStyle','-','Marker','+','MarkerIndices',tt);hold on;
    plot(t,error1(3,:),'Color','b','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,error(3,:)-error1(3,:),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
    % plot(t,U1(3,:),'Color','b','LineStyle','-.','Marker','+','MarkerIndices',tt);hold on;
    % plot(t,U2(3,:),'Color','g','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
    % plot(t,y_all(:,3),'Color','r','LineStyle','-','Marker','none','MarkerIndices',tt);hold on;
else
    figure,
    plot3(U(1,:),U(2,:),U(3,:),'g*');
    figure,
    plot3(U1(1,:),U1(2,:),U1(3,:),'g*');
end