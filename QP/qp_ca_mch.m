function u= qp_ca_mch(arg,B,plim,rlim,T,Wv,Wu,ud,imax,gam,only_plim)
  % �޸�qp_ca_slԴ�ļ��õ���ɾ�������㷨�ķ�֧��������matlabתc���Ժ�������飬��Ӧ��ϵΪ��
  % �����Ԫ������Ϊ�����һ�С��ڶ��С���������ֱ�����һ�С�
% Wrapper used in the QP control allocation Simulink block.
  
% Dimensions
  k=3; %B����
  m=9; %B����
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
u= wls_alloc_mch(B,v,umin,umax,Wv,Wu,ud,gam,uprev,W0,imax,m);