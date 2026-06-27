clear; clc;

% Sun et al. 2022 Table II equivalent physical scale.
m = 0.75;
J = diag([2.5, 2.1, 4.3])*1e-3;   % kg*m^2

CT = 8.5;                          % max thrust per normalized actuator [N]
cq = 2.37e-8;
ct = 1.51e-6;
kappa = cq/ct;                     % yaw moment / thrust [m]

% Iris Gazebo Classic geometry and actuator order, unchanged.
pos = [ ...
     0.13,  0.22, -0.023;
    -0.13, -0.20, -0.023;
     0.13, -0.22, -0.023;
    -0.13,  0.20, -0.023];
KM = [kappa; kappa; -kappa; -kappa];

axis = [0; 0; -1];

B_px4 = zeros(6,4);

for i = 1:4
    r = pos(i,:)';
    moment = CT * cross(r, axis) - CT * KM(i) * axis;
    force  = CT * axis;
    B_px4(:,i) = [moment; force];
end

[D_px4, B_norm_px4, mix_norm_px4, scale_px4, mix_raw_px4] = px4_normalize_B(B_px4, true);

B = [B_px4(6,:); B_px4(1:3,:)] % [Fz; Mx; My; Mz]
B_norm=[B_norm_px4(6,:); B_norm_px4(1:3,:)]
