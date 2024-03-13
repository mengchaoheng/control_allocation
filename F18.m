% F18 example
% clc;
% clear all;
% close all;
load f18data
[u,W,time,iter]=qp_sim(B,v,plim,rlim,T1);
figure,
plot(tn,u);
figure,
plot(tn,v,'k',tn,B*u);
