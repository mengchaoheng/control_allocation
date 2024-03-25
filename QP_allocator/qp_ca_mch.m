function u= qp_ca_mch(arg,B,plim,rlim,T,Wv,Wu,ud,imax,gam,only_plim)
  % 修改qp_ca_sl源文件得到，删除其他算法的分支。矩阵在matlab转c语言后会变成数组，对应关系为：
  % 数组的元素依次为矩阵第一列、第二列、。。。、直到最后一列。
% Wrapper used in the QP control allocation Simulink block.
  
% Dimensions
  % k=3; %B行数
  % m=4; %B列数
  [k,m] = size(B);
  % Extract nonconstant input arguments
  v = arg(1:k);
  uprev = arg(k+1:end);
  
  % Overall position limits
  if (only_plim)
    umin = plim(:,1);
    umax = plim(:,2);
  else
    umin = max(plim(:,1),uprev+rlim(:,1)*T);
    umax = min(plim(:,2),uprev+rlim(:,2)*T);
  end
% u0 = (umin+umax)/2;
% W0 = zeros(m,1);

% u0=uprev;
W0 = zeros(m,1);
infeasible1 = (uprev <= umin);
W0(infeasible1)=-1;
infeasible2 = (uprev >= umax);
W0(infeasible2)= 1;
%   u = wls_alloc(B,v,umin,umax,Wv,Wu,ud,gam,u0,W0,imax);
u= wls_alloc_gen(B,v,umin,umax,Wv,Wu,ud,gam,uprev,W0,imax,m);
