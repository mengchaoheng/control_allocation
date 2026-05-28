# 控制分配设计规律：SHC09 日志与 SHW09 VTOL 模式测试

本文总结 `test_shc09_px4_log_allocation.m` 和 `test_shw09_vtol_modes.m` 的主要结论，用来指导后续控制分配矩阵、限幅和失效模式设计。

对应测试文件：

- `test_shc09_px4_log_allocation.m`：PX4 日志、`B_norm/B_par`、合并分配与 split 分配对比。
- `test_shw09_vtol_modes.m`：VTOL 悬停/平飞模式、舵面置零/锁死/降维/完整体对比。

本文统一使用控制分配记号：

```text
y = B*u
```

其中 `y` 是控制目标，`B` 是 effectiveness matrix，`u` 是执行器输出。

## 1. PX4 INV 与 B_norm

PX4 的 `ControlAllocationPseudoInverse` 不是直接用原始物理 `B` 做 `pinv(B)*y`，而是：

1. 先对原始 effectiveness matrix 求伪逆：

   ```text
   M = B^+
   ```

2. 根据 `M` 的列计算 allocation scale：

   ```text
   D = diag(scale)
   ```

3. 对 mix 的列做归一化：

   ```text
   M_norm = M*D^-1
   ```

4. 实际输出：

   ```text
   u = M_norm*y
   ```

在活动控制轴行满秩、且 `D` 可逆时，这等价于先把 `B` 做行缩放：

```text
B_norm = D*B
```

然后用：

```text
u = pinv(B_norm)*y
```

注意这里的 `D` 来自 `B^+` 的列归一化，不是从 `B` 每一行直接求平均得到的行归一化。

对 SHC09 来说，完整 6 行控制轴为：

```text
[Mx My Mz Fx Fy Fz]
```

但 `Fx/Fy` 行全零，因此：

```text
rank(B_norm) = 4 < 6
```

所以用于 DP / PCA / direction-preserving 等优化算法时，应使用活动行：

```text
B_par = B_norm([Mx My Mz Fz], :)
```

这不是随意删维，而是删除 `0 = 0` 的退化约束。`inv/pinv` 对秩亏问题仍能给最小范数解，但很多 LP/PCA/DP 类方法依赖有效约束构造，完整零行会让问题退化。

## 2. B / B_norm 中哪些相对值重要

单位化以后，不能只看某个元素的绝对大小。更重要的是两个相对关系：

```text
同一行内各执行器的相对值，更直接决定这个控制轴由谁分担。
同一列内力/力矩的相对值，决定某个执行器带来的跨轴耦合和补偿。
```

例如 SHC09 中 motor0 同时产生 `Fz` 和 `Mz`。`Fz` 行决定 motor0 对推力的一对一增益；motor0 这一列里的 `Mz/Fz` 相对值决定它在产生推力时会附带多少偏航力矩，以及后续舵面要补偿多少 `Mz`。

## 3. 合并力和力矩分配

合并分配是直接把所有活动控制轴放进同一个问题：

```text
y_par = [Mx My Mz Fz]'
B_par = B_norm([Mx My Mz Fz], :)
u = allocator(B_par, y_par, umin, umax)
```

优点：

- 数学形式统一。
- 所有执行器耦合一次性进入优化。
- 如果算法目标正好就是全局目标，合并问题表达最完整。

风险：

- 对 DP / PCA / direction-preserving 类方法，合并问题里的冗余执行器和耦合列会带来多个等价最优顶点。
- LP 顶点选择可能改变 `u` 分担方式，即使 `B*u` tracking 很好，执行器曲线也可能和 PX4 inv 不同。
- 若某个执行器既产生力又产生力矩，例如 SHC09 motor0 的 `Fz + Mz`，合并分配可能为了方向保持或边界目标改变 motor0/servo 的分担。

SHC09 测试中，在修正 tie-break 后，`pca_dir_Bpar` 和 `inv_Bpar` 对日志窗口表现一致；但这依赖当前输入范围、约束状态和求解器 tie-break 稳定。合并问题仍然是更容易暴露 LP 顶点选择差异的形式。

## 4. 分开力和力矩分配

对 SHC09 这类结构，舵面不产生力，只有 motor0 产生 `Fz`。因此可以把问题写成：

```text
y = [tau; f]
u = [u_force; u_servo]

tau = B_torque_force*u_force + B_torque_servo*u_servo
f   = B_force_force*u_force + 0*u_servo
```

推荐流程：

1. 先计算力通道：

   ```text
   u_force_raw = pinv(B_force_force)*f
   u_force = clamp(u_force_raw, umin_force, umax_force)
   ```

   若是一对一通道，也可以直接除以增益。

2. 用限幅后的真实 `u_force` 计算它已经产生的力矩：

   ```text
   tau_from_force = B_torque_force*u_force
   ```

3. 从期望力矩里扣掉这部分：

   ```text
   tau_left = tau_cmd - tau_from_force
   ```

4. 用不产生力的执行器分配剩余力矩：

   ```text
   u_servo = allocator(B_torque_servo, tau_left, umin_servo, umax_servo)
   ```

5. 合成完整输出：

   ```text
   u = [u_force; u_servo]
   ```

关键点是：力通道限幅后，必须用限幅后的 `u_force` 去算补偿力矩。不要用未限幅的 raw 值。

这种做法适合 SHC09，因为：

- `B_force_servo = 0`；
- 力部分可独立先算；
- motor0 附带的 `Mz` 可以明确扣回 `tau_left`；
- 剩余力矩分配只在 3 维力矩空间里做，问题更干净。

## 5. 合并与分开的等价性

不要简单认为“同一个算法，先算一个再算一个”和“同时算”总是数学等价。

无约束伪逆时：

- 若 `B` 可通过行/列结构解耦成块对角，且目标函数范数权重也对应解耦，分开和合并可能得到相同结果。
- 若存在跨块耦合列，例如 `u_force` 同时产生 `Fz` 和 `Mz`，分开法等价于人为规定“力先满足，再用剩余执行器补偿力矩”。这通常是工程上想要的，但不一定等于整体 `pinv(B)` 的最小范数解。

有约束或限幅时：

- 一般不等价。
- 合并算法会在所有约束中寻找全局 LP/QP 解。
- 分开算法等价于加入优先级：力优先，力执行器先限幅，然后力矩补偿。
- 如果这是设计意图，分开更稳定、更可解释；如果希望所有轴统一折中，才应合并。

SHC09 当前日志窗口中，`inv_Bpar`、`split_inv`、`split_pca_dir`、`split_pca_dpscaled` 与 PX4 inv 日志输出非常接近，说明该构型和输入范围下 split 是干净可用的。更重要的是，split 避免了把 `Fz` 与力矩一起丢进 DP/PCA 时的冗余顶点问题。

## 6. 行满秩与 B_par

对优化类控制分配方法，尤其是 DP/PCA/LPCA：

- 输入矩阵应使用有效控制轴；
- 删除全零行；
- 尽量保证参与求解的 `B` 行满秩；
- 不要把“命令也为零的全零控制轴”留在 LP 约束里。

`inv/pinv` 可以处理秩亏，它会给一个最小范数解；但 DP/PCA 的 LP 构造、初始化问题和方向保持约束通常不喜欢这种退化行。

因此推荐命名：

```text
B_norm      = PX4 归一化后的完整坐标矩阵
B_par       = B_norm 删除全零行后的活动控制轴矩阵
B_torque    = B_par 中的力矩行
B_force     = B_par 中的力行
```

## 7. VTOL 悬停失效舵面的建模

SHW09 VTOL 测试比较了 5 类悬停/失效建模：

```text
1. B(:,7:8)=0 且 u7/u8 锁死为 0
2. 删除 u7/u8 列，降维为 6 执行器
3. B(:,7:8)=0，但 u7/u8 仍保持自由限幅
4. B(:,7:8) 保持非零，但 u7/u8 锁死为 0
5. 完整体：B(:,7:8) 非零，u7/u8 也自由
```

结论：

- 悬停阶段如果某些舵面实际不应参与分配，推荐同时做：

  ```text
  B(:, failed_or_disabled_idx) = 0
  umin(failed_or_disabled_idx) = 0
  umax(failed_or_disabled_idx) = 0
  ```

- 只把 `B` 列置零但不锁限幅，舵面变成纯零空间自由变量。部分算法或 restoring 可能仍移动这些舵，虽然 `B*u` 不变，但实际执行器会动，不适合真实悬停禁用。

- 只锁 `umin=umax=0`，但保留非零 `B`，是不推荐的反模式。求解器认为这些舵有控制效能，但又完全不可动，会导致 DP/PCA 类方法在其它执行器上产生奇怪顶点选择，restoring 也可能被锁死通道卡住。

- 降维删除列在数学上很干净，通常与“置零 B + 锁死限幅”的活动通道结果等价；但工程上会改变维度、索引、日志和下游接口，不一定方便。

因此 VTOL 模式切换推荐：

```text
悬停/舵面禁用:
    B(:, idx) = 0
    umin(idx) = 0
    umax(idx) = 0

平飞/舵面恢复:
    B(:, idx) = physical_effectiveness
    umin(idx) = physical_lower_limit
    umax(idx) = physical_upper_limit
```

也就是说，等舵面恢复效应后，再同时恢复对应 `B` 列和限幅。

## 8. restoring 的边界

`restoring_cpp` 只能沿 null-space 或可行方向整理执行器输出。它不是万能修复器。

典型失败场景：

- 求解器给出的 LP 顶点本身就是另一个等价顶点；
- 被锁死通道仍在 `B` 中有非零效能；
- restoring 的目标点想使用锁死通道；
- 锁死通道使恢复方向一开始就不可行。

SHW09 中 `mc-elevB-lock78` 就是反模式：`B(:,7:8)` 有舵效，但 `u7/u8` 被锁死。DP_LPCA 只最大化 `lambda`，没有“让 u 小一点/远离限幅”的二级目标；restoring 又会倾向完整 B 的伪逆目标，而那个目标想用 `u7/u8`。由于 `u7/u8` 被锁死，恢复过程会被卡住。

因此不要依赖 restoring 去弥补错误的 `B/limit` 建模。先把可用执行器和不可用执行器建模正确。

## 9. 推荐设计流程

1. 先明确坐标：

   ```text
   y = B*u
   ```

2. 如果要对齐 PX4 inv，先复现 PX4 scale：

   ```text
   M = px4_geninv(B)
   D = diag(scale_from_M)
   B_norm = D*B
   ```

3. 删除全零控制轴行，得到活动行：

   ```text
   B_par = B_norm(active_rows, :)
   y_par = y(active_rows)
   ```

4. 若力通道可以由一组执行器独立产生，且其它执行器不产生力，优先考虑 split：

   ```text
   先力，clamp，再扣除力执行器附带力矩，最后分配剩余力矩
   ```

5. 若执行器禁用或失效：

   ```text
   B 对应列置零
   umin = umax = 0
   ```

6. 如果降维实现成本低，可以删除禁用列；如果系统接口不想改维度，就保留完整维度并使用“置零 B + 锁死限幅”。

7. 不要只锁限幅而保留非零效能；不要只置零效能而让真实舵面自由乱动。

8. 只有在建模正确后，再用 restoring 优化执行器偏置或远离限幅。

## 10. 给工程实现的短结论

- PX4 inv 对比时，不要直接拿物理 `B` 做普通 `pinv(B)`；先构造 `B_norm = D*B`。
- DP/PCA/LPCA 类算法尽量喂行满秩的活动行 `B_par`。
- SHC09 这类“电机产生力并附带力矩、舵面只产生力矩”的构型，split 比合并更稳、更可解释。
- VTOL 悬停禁用舵面时，最佳工程做法是 `B` 列置零并同时锁死对应 `u`。
- 平飞恢复舵面时，必须同时恢复 `B` 列和限幅。
- restoring 不是失效建模的替代品；它只能在正确的可行集合内整理解。
