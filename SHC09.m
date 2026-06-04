clc    
% SHC09 six-effector model, same physical construction as alloc_cpp/test/main.cpp.
l1 = 0.267;
l2 = 0.066;
I_x = 0.0438;
I_y = 0.0436;
I_z = 0.005006;
d = 60*pi/180;
I = diag([I_x, I_y, I_z]);
Omega_max=1200;%
Omega_min=0;%
m=2.1;

R=0.114;
sigma_d=0.7268;

area_cs=0.00693;


% 先用悬停信息求等价参数
Omega_h=837; % 悬停油门0.71，转速1225  

motorConstant=m*9.8/(Omega_h*Omega_h); %  
C_T=motorConstant/(0.5*1.2041*pi*R^4); % 

% V_h_e=k_v * Omega;
V_h_e=sqrt(m*9.8/(sigma_d*1.2041*pi*R^2));%理论出口风速24，实测 
 
k_v=V_h_e/Omega_h; % parameter in propeller_wind_map
C_lcv=0.767;% px4: q=1/2 ρ V_h_e^2, control_joint_rad_to_cl=C_lcv
controlJointRadToCL_cs=C_lcv;

q_h_cs = 0.5 * 1.2041 * V_h_e^2;% 悬停 24 m/s   
dF_d_delta_h_cs=controlJointRadToCL_cs* q_h_cs * area_cs; %尾舵 理论   实测k_omega2force:  


k_h_tau=dF_d_delta_h_cs/V_h_e^2;  %悬停0.0032

k_cs = dF_d_delta_h_cs % k_cs=k_tau*V^2_e,  悬停 dF_d_delta_h_cs 或者巡航 dF_d_delta_fw_cs




%=============================6==================================
% Original geometry note:
% B = I\[-l1 -l1*cos(d) l1*cos(d) l1 l1*cos(d) -l1*cos(d);
%        0 -l1*sin(d) -l1*sin(d) 0 l1*sin(d) l1*sin(d);
%        l2 l2 l2 l2 l2 l2]*k_cs;
%      1(X)
%  2       6
% -------------->(Y)
%  3       5
%      4
%
% The active sign convention below matches alloc_cpp/test/main.cpp.
B = I \ [-l1, -l1*cos(d),  l1*cos(d),  l1,  l1*cos(d), -l1*cos(d);
          0,   l1*sin(d),  l1*sin(d),  0,  -l1*sin(d), -l1*sin(d);
          l2,  l2,         l2,         l2,  l2,         l2] * k_cs
K = I\diag([l1 l1 l2])*k_cs;
P = [-1     -cos(d)     cos(d)    1    cos(d)  -cos(d);
      0,    sin(d),     sin(d),   0,  -sin(d), -sin(d);
      1     1           1         1    1        1]%B = K*P
% B_inv_pid=[-1 0 1;-1 1 1;1 1 1;1 0 1;1 -1 1;-1 -1 1]
% B_pid=pinv(B_inv_pid)=[-0.1667   -0.1667    0.1667    0.1667    0.1667   -0.1667;
%                         0         0.2500    0.2500    0        -0.2500   -0.2500;
%                         0.1667    0.1667    0.1667    0.1667    0.1667    0.1667];
ulim = 40*pi/180;



% 对应 PX4 内置 PseudoInverse/SequentialDesaturation 的实际单位化过程：
% updateControlAllocationMatrixScale() + normalizeControlAllocationMatrix()。
[D, B_norm, mix_norm, scale, mix_raw] = px4_normalize_B(B, true)


% PX4:
% pxh> control_allocator status
% INFO  [control_allocator] Running
% INFO  [control_allocator] Method: Auto
% INFO  [control_allocator] Effectiveness Source: Custom
% INFO  [control_allocator]   Effectiveness =
%   | 0      | 1      | 2      | 3      | 4      | 5      | 6      
% Mx| 0       -0.50100 -0.25050  0.25050  0.50100  0.25050 -0.25050 
% My| 0        0        0.43390  0.43390  0       -0.43390 -0.43390 
% Mz| 0        0.20700  0.20700  0.20700  0.20700  0.20700  0.20700 
% Fx| 0        0        0        0        0        0        0       
% Fy| 0        0        0        0        0        0        0       
% Fz|-6.50000  0        0        0        0        0        0       
% INFO  [control_allocator]   minimum =
%   | 0      | 1      | 2      | 3      | 4      | 5      | 6      
%   | 0       -1.00000 -1.00000 -1.00000 -1.00000 -1.00000 -1.00000 
% INFO  [control_allocator]   maximum =
%   | 0      | 1      | 2      | 3      | 4      | 5      | 6      
%   | 1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000 
% INFO  [control_allocator]   Configured actuators: 7
% control_allocator: cycle: 6594 events, 0us elapsed, 0.00us avg, min 0us max 0us 0.000us rms
% pxh> 
% 舵机（力矩）部分
B_px4=[-0.50100 -0.25050  0.25050  0.50100  0.25050 -0.25050;
        0        0.43390  0.43390  0       -0.43390 -0.43390;
        0.20700  0.20700  0.20700  0.20700  0.20700  0.20700]
% 旧 B_px4 也按 PX4 内置 RPY normalization 算尺度，便于和当前物理 B 对齐。
[D_px4, B_norm_px4, mix_norm_px4, scale_px4, mix_raw_px4] = px4_normalize_B(B_px4, true)
% U = (1/umax) * E单位矩阵 
% D * B = D_px4 * B_px4, 
% INDI 控制器增益
% k = D^-1 * k_px4 / U，即
k_px4 =  [0.4;0.4;0.3];

k_indi= inv(D) * k_px4 * ulim
% 内置分配非单位化的情形，还需要* inv(D_px4)
k_indi = inv(D) * inv(D_px4) * k_px4 * ulim
