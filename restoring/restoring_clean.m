function u_rest = restoring_clean(B, u, uMin, uMax)

tol = 1e-5;

[~,m] = size(B);

% 1. 当前控制效果对应的最小范数解
u_pseudo = B' * ((B*B') \ (B*u));

% 2. 当前 u 的零空间分量
u_n = u - u_pseudo;
n2 = u_n' * u_n;

if n2 < tol^2
    u_rest = u;
    return;
end

% 3. 恢复方向：朝 -u_n
d = -u_n;

% 4. 最优步长：若无限幅，K=1 直接到 u_pseudo
K_opt = 1;

% 5. 限幅允许的最大步长
K_max = Inf;
for i = 1:m
    if abs(d(i)) < tol
        continue;
    end

    if d(i) > 0
        tmp = (uMax(i) - u(i)) / d(i);
    else
        tmp = (uMin(i) - u(i)) / d(i);
    end

    K_max = min(K_max, tmp);
end

K = min(K_opt, K_max);
K = max(0, K);

u_rest = u + K*d;

end