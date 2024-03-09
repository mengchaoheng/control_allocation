function u = dyn_alloc_m(B,v,u,p_limits,v_limits,T,Wv,W1,W2,S,imax,gam,only_plim)
  % Formulate as a regular QP problem: min ||Wu(u-ud)||
  W1sq = W1^2;
  W2sq = W2^2;
  % 非对角阵用sqrtm
%   Wu = real(sqrtm(W1sq+W2sq));
 % 对角阵用sqrt
  Wu = sqrt(W1sq+W2sq);
  us = S*v;
  ud = (W1sq+W2sq)\(W1sq*us+W2sq*u);
  % Use WLS -- fast and robust.
  u = wls_alloc_m(B,v,u,p_limits,v_limits,T,Wv,Wu,ud,imax,gam,only_plim);