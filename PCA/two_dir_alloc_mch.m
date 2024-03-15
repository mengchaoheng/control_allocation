 function [u]=two_dir_alloc_mch(v_T, v_D,p_limits,v_limits,u)
% (c) mengchaoheng
% Last edited 2019-11
v=v_T+v_D;
umin=[1;1;1;1]*(-p_limits)*pi/180;
umax=[1;1;1;1]*p_limits*pi/180;
% uv = allocator_dir_simplex_4(v,umin,umax); % wls_alloc_mch(v,u);% 先计算合力矩所需舵量
uv =wls_ca_4df_pv_limit(v, u,p_limits, v_limits);
if(check(uv,p_limits)) % 若舵量可以满足则直接满足
    u=uv;
else  % 否则再计算扰动所需
%     uv1 = allocator_dir_simplex_4(v_T,umin,umax); % wls_alloc_mch(v1,u);
    uv1 =wls_ca_4df_pv_limit(v_T,u, p_limits, v_limits);
    if(check(uv1,p_limits))  % 若扰动可满足，合力矩不能满足，则进行两次分配
        umin1=umin-uv1;
        umax1=umax-uv1;
        uv2 = allocator_dir_simplex_4(v_D,umin1,umax1);
%         uv2 =wls_ca_4df_pv_limit(v_D,uv1, umin1, umax1);
        u=uv1+uv2;
    else  % 扰动也不能满足，则直接按照合力矩进行分配
        u=uv;
    end
end
