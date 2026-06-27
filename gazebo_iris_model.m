clear; clc

% ===== 控制分配输出 =====
% u_alloc 是 PX4 control allocator 输出，normalized thrust
u_alloc = [0.5; 0.5; 0.5; 0.5];
u_alloc = min(max(u_alloc, 0), 1);

% ===== PX4 THR_MDL_FAC: normalized thrust -> motor signal =====
% PX4 模型: u_alloc = a*x^2 + (1-a)*x
% x 是送到仿真输出层的 normalized motor signal
a = 1.0;  % THR_MDL_FAC

a = min(max(a, 0), 1);

if abs(a) < 1e-9
    x = u_alloc;
else
    x = (-(1-a) + sqrt((1-a)^2 + 4*a*u_alloc)) / (2*a);
end

x = min(max(x, 0), 1);

% ===== Gazebo Classic Iris: motor signal -> commanded rotor speed =====
% PX4/SITL -> Gazebo mavlink_interface
% armed:
% omega = (x + input_offset) * input_scaling + zero_position_armed
omega = 1000*x + 0;
omega = min(max(omega, 0), 1000);

% ===== SDF motor model: rotor speed -> thrust =====
% SDF: F = motorConstant * omega^2
motorConstant = 8.5e-06;
F = motorConstant * omega.^2;

% ===== SDF physical B: motor thrust -> body wrench =====
% y = [Mx; My; Mz; Fx; Fy; Fz]
% 输入是每个电机的实际推力 F_i，单位 N

pos = [ ...
     0.13,  0.22, -0.023;
    -0.13, -0.20, -0.023;
     0.13, -0.22, -0.023;
    -0.13,  0.20, -0.023];

spin = [1; 1; -1; -1];   % ccw=+1, cw=-1
momentConstant = 0.0156953642384106;

B_force = zeros(6,4);

for i = 1:4
    r = pos(i,:)';
    force_per_N = [0; 0; -1];
    moment_per_N = cross(r, force_per_N) + [0; 0; spin(i)*momentConstant];

    B_force(:,i) = [moment_per_N; force_per_N];
end

% ===== 实际物理力和力矩 =====
y = B_force * F;

% Physical SDF model:
%   omega_cmd = input_scaling*x + zero_position_armed
%   F = motorConstant * omega_cmd^2
%
% For allocation we choose/interpret the actuator mapping so that:
%   F ~= Fmax * u_alloc
%
% Then:
%   y = B_force * F
%     = B_force * Fmax * u_alloc
%     = B_px4 * u_alloc

F_max=motorConstant * 1000.^2;
B_px4= B_force * F_max 

% B_px4 =
% 
%    -1.2848    1.1680    1.2848   -1.1680
%     0.7592   -0.7592    0.7592   -0.7592
%     0.2920    0.2920   -0.2920   -0.2920
%          0         0         0         0
%          0         0         0         0
%    -5.8400   -5.8400   -5.8400   -5.8400