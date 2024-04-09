function u=prioritized_control_allocator(m1,m2,B,umin,umax,itlim)
    [u_all, errout_all, lambda_all] = DPscaled_LPCA(m1+m2,B,umin,umax,itlim);
    % the problem above always have solution
    if(lambda_all<1)
        % disp('m1+m2不可达, 计算m1可达性');
        [u_m1, errout_m1, lambda_m1] = DPscaled_LPCA(m1,B,umin,umax,itlim);
        if(lambda_m1<1)
            % disp('m1不可达, 重新构造问题');
            % disp('使用m1+m2的分配结果, 从m1+m2方向收缩');
            % disp('最终力矩：');
            % B*u_all
            u=u_all;
        else 
            % disp('m1可达，重新计算');
            [u_m2, errout_m2, lambda_m2] = DPscaled_LPCA(m2,B,umin-u_m1,umax-u_m1,itlim);
            % disp('最终分配结果：');
            u=u_m1+u_m2;
            % disp('最终力矩：');
            % B*(u_m1+u_m2)
        end
    
    else
        u=u_all;
        % disp('m1+m2可达, 从m1+m2方向收缩');
        % disp('最终分配结果：');
        % u_all
        % disp('最终力矩：');
        % B*u_all
        % disp('如果从m2方向收缩，会出现不符合要求的值，特别是m1不可达时');
    end
end