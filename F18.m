clear all;
close all;
addpath(genpath(pwd))
folder ='some_modified_function'; 
rmpath(folder) % remove old version
folder ='s-function_used_in PlanD'; 
rmpath(folder) % remove old version

load f18data
[u,W,time,iter]=qp_sim(B,v,plim,rlim,T1,'alg','sls');
figure,
plot(tn,u);
figure,
plot(tn,v,'k',tn,B*u);
