function [isOuted, out_mask] = isOut(u, lower_bound, upper_bound)
% 判断向量 u 是否超过边界 
% 输入：
%   u: m 维输入向量
%   lower_bound: 下界标量或与 u 同维度的向量
%   upper_bound: 上界标量或与 u 同维度的向量
% 输出：
%   is_saturated: 布尔值。若 u 中至少一个元素超过边界则为 true，否则 false
%   saturated_mask: 布尔向量。标记 u 中每个元素是否超过边界（1=超过边界，0=容许）

% 输入参数校验
if nargin < 3
    error('需提供下界 lower_bound 和上界 upper_bound');
end
if any(lower_bound > upper_bound)
    error('下界不能超过上界');
end

% 扩展标量边界为向量（若输入为标量）
if isscalar(lower_bound)
    lower_bound = lower_bound * ones(size(u));
end
if isscalar(upper_bound)
    upper_bound = upper_bound * ones(size(u));
end
% e表示一点点容忍误差。
e=100*eps;
% 检测饱和元素：u < lower_bound - e 或 u > upper_bound + e 
below_lower = (u < lower_bound-e);
above_upper = (u > upper_bound+e);
out_mask = below_lower | above_upper;

% 判断整体是否超过边界
isOuted = any(out_mask(:));
end