clear; clc;


I_x = 0.01149;
I_y = 0.01153;
I_z = 0.00487;
l1 = 0.167;
l2 = 0.069;
I = diag([I_x, I_y, I_z]);
Omega_max=1750;%
Omega_min=0;%
m=1.39;

R=0.114;
sigma_d=0.7;

area_cs=0.00693;


% 先用悬停信息求等价参数
Omega_h=1225; % 悬停油门0.71，转速1225  

motorConstant=m*9.8/(Omega_h*Omega_h); %  
C_T=motorConstant/(0.5*1.2041*pi*R^4); % 

% V_h_e=k_v * Omega;
V_h_e=sqrt(m*9.8/(sigma_d*1.2041*pi*R^2));%理论出口风速19.8958，实测 
 
k_v=V_h_e/Omega_h; % parameter in propeller_wind_map
C_lcv=1.7266;% px4: q=1/2 ρ V_h_e^2, control_joint_rad_to_cl=C_lcv
controlJointRadToCL_cs=C_lcv;

q_h_cs = 0.5 * 1.2041 * V_h_e^2;% 悬停 24.1 m/s   
dF_d_delta_h_cs=controlJointRadToCL_cs* q_h_cs * area_cs; %尾舵 理论4.2055  实测k_omega2force: 4.203


k_h_tau=dF_d_delta_h_cs/V_h_e^2;  %悬停0.0072 

k_cs = dF_d_delta_h_cs % k_cs=k_tau*V^2_e,  悬停 dF_d_delta_h_cs 或者巡航 dF_d_delta_fw_cs



B = I \ [-l1     0       l1     0;
          0      -l1     0      l1;
          l2     l2      l2     l2] * k_cs
ulim = 20*pi/180;


% 对应 PX4 内置 PseudoInverse/SequentialDesaturation 的实际单位化过程：
% updateControlAllocationMatrixScale() + normalizeControlAllocationMatrix()。
[D, B_norm, mix_norm, scale, mix_raw] = px4_normalize_B(B, true)
 
% PX4 control_allocator 打印出来的 DF4 悬停分配矩阵。
% pxh> control_allocator status
% INFO  [control_allocator] Running
% INFO  [control_allocator] Method: Auto
% INFO  [control_allocator] Effectiveness Source: Custom
% INFO  [control_allocator] B Unit: User-defined
% INFO  [control_allocator] Allocation Unit: Normalized
% INFO  [control_allocator]   Effectiveness =
%   | 0      | 1      | 2      | 3      | 4      
% Mx| 0       -0.50000  0        0.50000  0       
% My| 0        0       -0.50000  0        0.50000 
% Mz| 0.32500  0.25000  0.25000  0.25000  0.25000 
% Fx| 0        0        0        0        0       
% Fy| 0        0        0        0        0       
% Fz|-6.50000  0        0        0        0       
% INFO  [control_allocator]   minimum =
%   | 0      | 1      | 2      | 3      | 4      
%   | 0       -1.00000 -1.00000 -1.00000 -1.00000 
% INFO  [control_allocator]   maximum =
%   | 0      | 1      | 2      | 3      | 4      
%   | 1.00000  1.00000  1.00000  1.00000  1.00000 
% INFO  [control_allocator]   Configured actuators: 5
% control_allocator: cycle: 3462 events, 0us elapsed, 0.00us avg, min 0us max 0us 0.000us rms
% pxh> 
% 控制量顺序: y = [Mx; My; Mz]
% 执行器顺序: u = [servo0 servo1 servo2 servo3]'
B_px4 = [ -0.5  0    0.5  0 ;
           0   -0.5  0    0.5;
           0.25 0.25 0.25 0.25];
% 旧 B_px4 也按 PX4 内置 RPY normalization 算尺度，便于和当前物理 B 对齐。
[D_px4, B_norm_px4, mix_norm_px4, scale_px4, mix_raw_px4] = px4_normalize_B(B_px4, true)


 

% D * B = D_px4 * B_px4,  
% U = (1/umax) * E单位矩阵 
% 分配模型 v=Bu -> D * v * U = D * B * U * u_norm = B_norm * u_norm 

% px4: v_px4 = k_px4 * e = D_px4 * B_px4 * u_px4 = B_norm * u_px4 (px4控制器输出v_px4已经单位化，操纵面u_px4也默认单位化，实际用于分配的是B_norm=D_px4 * B_px4)
% INDI: v=k_indi * e,  D * v * U = D * B * U * u_indi = B_norm * u_indi  (INDI单位化需要两边同乘D,还有u_indi的单位化产生的U)
% 要使u_indi接近u_px4
% D * k_indi * e * U = k_px4 * e
% 因此控制器增益
% k = D^-1 * k_px4 / U，即
k_px4 =  [0.45;0.45;0.3];
k_indi = inv(D) * k_px4 * ulim  
% 内置分配非单位化的情形，还需要* inv(D_px4)
k_indi = inv(D) * inv(D_px4) * k_px4 * ulim
