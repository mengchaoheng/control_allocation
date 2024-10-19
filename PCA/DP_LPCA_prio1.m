function [u, errout,lambda] = DP_LPCA_prio1(m_higher,m_lower,B,uMin,uMax,itlim)
[k,~] = size(B);

[u, errout, lambda] = DP_LPCA_copy(m_higher,m_lower,B,uMin,uMax,itlim);
if(errout<0)
    % 构造新问题
    % disp('errout==-2无解，即m1+m2不可达且m1不可达');% 不可通过单独收缩m2实现
    % disp('or errout==-1 m2=0');% 
    [u, errout, lambda] = DP_LPCA_copy(zeros(k,1),m_higher,B,uMin,uMax,itlim);
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