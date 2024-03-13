function u = dyn_ca_4df(v,u)
%============for df which have 4 control surface================
B=[-0.5   0       0.5   0;
    0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
p_limits=20*pi/180;% rad
k=3;
m=4;
% [k,m] = size(B);
W2=1*eye(m);% maybe every value is different 
W1=1*eye(m);% 
plim(:,1)=ones(m,1)*(-p_limits);
plim(:,2)=ones(m,1)*p_limits;
rlim=[]; % rad/s
S=pinv(B);
%===============================================================
uprev=u;
Wv=eye(k);
% Formulate as a regular QP problem: min ||Wu(u-ud)||
W1sq = W1^2;
W2sq = W2^2;
% diagonal matrix use sqrt
Wu = sqrt(W1sq+W2sq);
% else sqrtm
%   Wu = real(sqrtm(W1sq+W2sq));
us = S*v;
ud = (W1sq+W2sq)\(W1sq*us+W2sq*uprev);
% Overall position limits
if isempty(rlim)
    umin = plim(:,1);
    umax = plim(:,2);
else
    umin = max(plim(:,1),uprev+rlim(:,1)*T);
    umax = min(plim(:,2),uprev+rlim(:,2)*T);
end
gam=1e6;
u0 = (umin+umax)/2;
W0 = zeros(m,1);
% number of iter
imax=100;
% Use WLS -- fast and robust.
u= wls_alloc_gen(B,v,umin,umax,Wv,Wu,ud,gam,u0,W0,imax,m);