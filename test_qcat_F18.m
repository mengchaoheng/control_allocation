clear all;
close all;
addpath(genpath(pwd))


load f18data
[u,W,time,iter]=qp_sim(B,v,plim,rlim,T1,'alg','sls');
figure,
plot(tn,u);
figure,
plot(tn,v,'k',tn,B*u);
