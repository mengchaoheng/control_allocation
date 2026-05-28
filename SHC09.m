clc    
% SHC09 six-effector model, same physical construction as alloc_cpp/test/main.cpp.
l1 = 0.292166;
l2 = 0.073699;
k_omega2force = 1.93;
I_x = 0.0438;
I_y = 0.0436;
I_z = 0.005006;
d = 60*pi/180;
I = diag([I_x, I_y, I_z]);

%=============================6==================================
% Original geometry note:
% B = I\[-l1 -l1*cos(d) l1*cos(d) l1 l1*cos(d) -l1*cos(d);
%        0 -l1*sin(d) -l1*sin(d) 0 l1*sin(d) l1*sin(d);
%        l2 l2 l2 l2 l2 l2]*k_v;
%      1(X)
%  2       6
% -------------->(Y)
%  3       5
%      4
%
% The active sign convention below matches alloc_cpp/test/main.cpp.
B = I \ [-l1, -l1*cos(d),  l1*cos(d),  l1,  l1*cos(d), -l1*cos(d);
          0,   l1*sin(d),  l1*sin(d),  0,  -l1*sin(d), -l1*sin(d);
          l2,  l2,         l2,         l2,  l2,         l2] * k_omega2force
K = I\diag([l1 l1 l2])*k_omega2force;
P = [-1     -cos(d)     cos(d)    1    cos(d)  -cos(d);
      0,    sin(d),     sin(d),   0,  -sin(d), -sin(d);
      1     1           1         1    1        1]%B = K*P
% B_inv_pid=[-1 0 1;-1 1 1;1 1 1;1 0 1;1 -1 1;-1 -1 1]
% B_pid=pinv(B_inv_pid)=[-0.1667   -0.1667    0.1667    0.1667    0.1667   -0.1667;
%                         0         0.2500    0.2500    0        -0.2500   -0.2500;
%                         0.1667    0.1667    0.1667    0.1667    0.1667    0.1667];
ulim = 40*pi/180;
umin = ones(6, 1) * -ulim;
umax = ones(6, 1) * ulim;


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

D = diag(scale)
B_norm = D * B


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
mix_raw_px4 = pinv(B_px4);          % MATLAB 近似 PX4 geninv；这个 B 下通常足够接近
scale_px4 = ones(size(B_px4,1),1);  % SHW09 这里按各控制轴列计算

for j = 1:size(B,1)
    col = mix_raw_px4(:,j);
    nz = abs(col) > eps_col;

    if any(nz)
        scale_px4(j) = mean(abs(col(nz)));
    else
        scale_px4(j) = 1;
    end
end

D_px4 = diag(scale_px4)
B_norm_px4 = D_px4 * B_px4

% D * B = D_px4 * B_px4, B = inv(D) * D_px4 * B_px4
% INDI 控制器增益
% k = U * D^-1 * k_norm，即
k_indi= inv(D) * D_px4 * [0.4;0.4;0.3] / ulim


 