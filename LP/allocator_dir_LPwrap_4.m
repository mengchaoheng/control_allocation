function [u,z,iters] = allocator_dir_LPwrap_4(B, v, umin,umax)
%%  come from [u] = LPwrap(IN_MAT) and using sigle data, define k, m for df4
% [k,m]=size(B);

itlim=5e2;

[u,itlim,~,z] = DPscaled_LPCA_C(v,B,umin,umax,itlim);
iters=5e2-itlim;
end