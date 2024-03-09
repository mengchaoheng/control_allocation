function [u]=two_dir_alloc_mch(v_T, v_D,p_limits,v_limits,u)
% (c) mengchaoheng
% Last edited 2019-11
v=v_T+v_D;
umin=[1;1;1;1]*(-p_limits)*pi/180;
umax=[1;1;1;1]*p_limits*pi/180;
% uv = dir_alloc_mch(v,umin,umax); % wls_alloc_mch(v,u);% �ȼ���������������
uv =wls_alloc_mch(v, u,p_limits, v_limits);
if(check(uv,p_limits)) % ����������������ֱ������
    u=uv;
else  % �����ټ����Ŷ�����
%     uv1 = dir_alloc_mch(v_T,umin,umax); % wls_alloc_mch(v1,u);
    uv1 =wls_alloc_mch(v_T,u, p_limits, v_limits);
    if(check(uv1,p_limits))  % ���Ŷ������㣬�����ز������㣬��������η���
        umin1=umin-uv1;
        umax1=umax-uv1;
        uv2 = dir_alloc_mch(v_D,umin1,umax1);
%         uv2 =wls_alloc_mch(v_D,uv1, umin1, umax1);
        u=uv1+uv2;
    else  % �Ŷ�Ҳ�������㣬��ֱ�Ӱ��պ����ؽ��з���
        u=uv;
    end
end
