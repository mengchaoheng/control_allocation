clear; clc;

I_x = 0.050636;
I_y = 0.042954;
I_z = 0.012668;
L_1 = 0.2399;
L_2 = 0.0664;
L_3 = 0.3145;
Omega_max=1555;%
Omega_min=466;%
m=2.05;
C_T0=0.0307;
R=0.114;
sigma_d=0.7;
area_wing=0.115;
area_cs=0.00693;

% 先用悬停信息求等价参数
Omega_h=0.7*(Omega_max-Omega_min)+Omega_min; % 悬停油门0.7，转速1228.3   
% m=(Omega_h*Omega_h)*(C_T0*0.5*1.2041*pi*R^4)/9.8 % 1.5 
% V_omega_h=sqrt(m*9.8/(C_T0*0.5*1.2041*pi*R^4)); % 1431

motorConstant=m*9.8/(Omega_h*Omega_h); % (C_T0*0.5*1.2041*pi*R^4) = k_T
C_T=motorConstant/(0.5*1.2041*pi*R^4); % 0.0417 说明实际的系数C_T区别于CT0

zero_position_armed=Omega_min;
input_scaling=Omega_max-Omega_min;

% V_h_e=k_v * Omega;
V_h_e=sqrt(m*9.8/(sigma_d*1.2041*pi*R^2));%理论出口风速24.1618，实测[24.1547]
k_v=V_h_e/Omega_h; % parameter in propeller_wind_map
C_lcv=1.7266;% px4: q=1/2 ρ V_h_e^2, control_joint_rad_to_cl=C_lcv
controlJointRadToCL_cs=C_lcv;

q_h_cs = 0.5 * 1.2041 * V_h_e^2;% 悬停 24.1 m/s   
dF_d_delta_h_cs=controlJointRadToCL_cs* q_h_cs * area_cs; %尾舵 理论4.2055  实测k_omega2force: 4.203

% 巡航参数
controlJointRadToCL_wing=0.3;%估计，论文无直接出处
V_a=16.4911; % 实测空速
V_fw_e=32.7686; % 实测出口风速
Omega_fw=0.34*(Omega_max-Omega_min)+Omega_min;% 实测巡航油门0.34
 
V_fw_e_omega=k_v*Omega_fw; %理论滑流 16.4500 实测16.2775=V_fw_e-V_a

q_fw_cs = 0.5 * 1.2041 * (V_fw_e)^2;%   
q_wing = 0.5 * 1.2041 * V_a^2;%  

dF_d_delta_fw_cs=controlJointRadToCL_cs* q_fw_cs * area_cs;%尾舵 理论7.7352  实测7.73524 
dF_d_delta_wing=controlJointRadToCL_wing*q_wing * area_wing;%机翼 理论5.6487  实测5.64871

k_e=dF_d_delta_wing/V_a^2; %k_e=0.0208
% k_tau = k_h_tau = k_fw_tau 系数相同，不同在于出口风速
k_h_tau=dF_d_delta_h_cs/V_h_e^2;  %悬停0.0072 
k_fw_tau=dF_d_delta_fw_cs/V_fw_e^2; % 巡航 0.0072 

k_wing = dF_d_delta_wing;  % k_wing= k_e * V^2_a
k_cs = dF_d_delta_h_cs; % k_cs=k_tau*V^2_e,  悬停 dF_d_delta_h_cs 或者巡航 dF_d_delta_fw_cs

d = 60*pi/180;
ulim = 0.6981;  

B = zeros(3, 8);
B(1,1:6) = [-L_1, -cos(d)*L_1, cos(d)*L_1, L_1, cos(d)*L_1, -cos(d)*L_1] * k_cs / I_x;
B(2,1:6) = [0, sin(d)*L_1, sin(d)*L_1, 0, -sin(d)*L_1, -sin(d)*L_1] * k_cs / I_y;
B(3,1:6) = ones(1,6) * L_2 * k_cs / I_z;

% 若悬停时 u7/u8 不参与，保持 B(:,7:8)=0 即可。
% 若固定翼要恢复 elevon，把 B(3,7:8) 填成对应非零值。
B(3,7) =  k_wing * L_3 / I_z; % 论文中 k_wing= k_e * V^2_a, 这里隐去风速
B(3,8) = -k_wing * L_3 / I_z;

mix_raw = pinv(B);          % MATLAB 近似 PX4 geninv；这个 B 下通常足够接近
scale = ones(size(B,1),1);  % SHW09 这里按各控制轴列计算

eps_col = 1e-6;
for j = 1:size(B,1)
    col = mix_raw(:,j);
    nz = abs(col) > eps_col;

    if any(nz)
        scale(j) = mean(abs(col(nz)));
    else
        scale(j) = 1;
    end
end

D = diag(scale);
B_norm = D * B;

% PX4 control_allocator 打印出来的 SHW09 悬停分配矩阵。
% 控制量顺序: y = [Mx; My; Mz]
% 执行器顺序: u = [servo0 servo1 servo2 servo3 servo4 servo5 elevon0 elevon1]'
B_px4 = [
   -0.50100  -0.25050   0.25050   0.50100   0.25050  -0.25050   0        0
    0         0.43390   0.43390   0        -0.43390  -0.43390   0        0
    0.20700   0.20700   0.20700   0.20700   0.20700   0.20700   0.50000 -0.50000
];

mix_raw_px4 = pinv(B_px4);
scale_px4 = ones(size(B_px4,1),1);

for j = 1:size(B_px4,1)
    col = mix_raw_px4(:,j);
    nz = abs(col) > eps_col;

    if any(nz)
        scale_px4(j) = mean(abs(col(nz)));
    else
        scale_px4(j) = 1;
    end
end

D_px4 = diag(scale_px4);
B_norm_px4 = D_px4 * B_px4;

disp('=== B ===');
disp('B =');
disp(B);

disp('=== B normalization ===');
disp('scale =');
disp(scale.');
  
disp('B_norm = D * B =');
disp(B_norm);


% 

disp('=== B_px4 ===');
disp('B_px4 =');
disp(B_px4);

disp('=== B_px4 normalization ===');
disp('scale_px4 =');
disp(scale_px4.');

disp('B_norm_px4 = D_px4 * B_px4 =');
disp(B_norm_px4);
 
k_indi= inv(D) * D_px4 * [0.5;0.5;0.4]/ulim