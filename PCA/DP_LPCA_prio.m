function [u, errout, lambda] = DP_LPCA_prio(m_higher, m_lower, B, uMin, uMax, itlim, varargin)
% Prioritized DP_LPCA.
%
% This keeps the original DP_LPCA_prio name and behavior.  The lower level
% DP_LPCA_copy uses the original book/MATLAB simplxuprevsol.m backend.
% Extra varargin inputs are accepted only for legacy scripts and are ignored.

[k, ~] = size(B);

[u, errout, lambda] = DP_LPCA_copy(m_higher, m_lower, B, uMin, uMax, itlim);
if errout < 0
    % 构造新问题
    % errout==-2 无解，即 m1+m2 不可达且 m1 不可达；不可通过单独收缩 m2 实现。
    % errout==-1 可由 m_lower=0 触发，因为 DP_LPCA_copy 检查的是 y_d。
    [u, errout, lambda] = DP_LPCA_copy(zeros(k, 1), m_higher, B, uMin, uMax, itlim); % 第二个参数取 m_higher 或者 m_higher+m_lower
% else
    % disp('有解，m1可达或者m1+m2可达');
    % disp('原问题解产生的力矩总是满足要求');
    % m_real=B*u
    % disp('解统一表达式是m1+lambda*m2');
    % m1+lambda*m2
    % if(lambda<1)
    %     disp('m1可达（无论m1+m2取何值总能有解），但lamb<1表示m1+m2不可达');
    %     disp('解是m1+lambda*m2:');
    %     m1+lambda*m2
    % else
    %     disp('m1+m2可达（无论m1），lamb=1');
    %     disp('解是m1+m2:');
    %     m1+m2
    % end
end

end
