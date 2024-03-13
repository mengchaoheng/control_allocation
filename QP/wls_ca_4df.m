function [u] = wls_ca_4df(v, u)
%============for df which have 4 control surface================
B=[-0.5   0       0.5   0;
    0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
p_limits=20*pi/180;% rad
k=3;
m=4;
% [k,m] = size(B);

plim(:,1)=ones(m,1)*(-p_limits);
plim(:,2)=ones(m,1)*p_limits;
rlim=[]; % rad/s
uprev=u;
%%
if isempty(rlim)
    umin = plim(:,1);
    umax = plim(:,2);
else
    umin = max(plim(:,1),uprev+rlim(:,1)*T);
    umax = min(plim(:,2),uprev+rlim(:,2)*T);
end
%%
Wv=eye(k);
Wu=eye(m);
ud=[0;0;0;0];% desire u
gam=1e6;
% u0 = mean([umin umax]')'; % init u, can be use uprev or 0
u0 = (umin+umax)/2;
W0 = zeros(m,1);
%--------
% if use u0=uprev; then the following code is need
% W0 = zeros(m,1);
% infeasible1 = (uprev <= umin);
% W0(infeasible1)=-1;
% infeasible2 = (uprev >= umax);
% W0(infeasible2)= 1;
%--------
% number of iter
imax=100;

u = wls_alloc_gen(B,v,umin,umax,Wv,Wu,ud,gam,u0,W0,imax,m);
end