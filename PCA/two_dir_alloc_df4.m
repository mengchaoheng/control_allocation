function [u]=two_dir_alloc_df4(v_T, v_D,umin,umax)
% (c) mengchaoheng
% Last edited 2019-11
v=v_T+v_D;
plim=[umin umax];
uv = allocator_dir_simplex_4(v,umin,umax); % wls_alloc_mch(v,u);% 先计算合力矩所需舵量
if(~check_if_saturation(uv,plim)) % 若舵量可以满足则直接满足
    u=uv;
else  % 否则再计算扰动所需
    uv1 = allocator_dir_simplex_4(v_T,umin,umax); % wls_alloc_mch(v1,u);
    if(~check_if_saturation(uv1,plim))  % 若扰动可满足，合力矩不能满足，则进行两次分配
        umin1=umin-uv1;
        umax1=umax-uv1;
        uv2 = allocator_dir_simplex_4(v_D,umin1,umax1);
        u=uv1+uv2;
    else  % 扰动也不能满足，则直接按照合力矩进行分配
        u=uv;
    end
end
