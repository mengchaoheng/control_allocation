function [D, B_norm, mix_norm, scale, mix_raw] = px4_normalize_B(B, normalize_rpy)
%PX4_NORMALIZE_B Reproduce PX4 PseudoInverse allocation normalization.
%   PX4 first computes mix_raw = pinv(B). If metric allocation is disabled,
%   updateControlAllocationMatrixScale() computes scale from mix_raw, then
%   normalizeControlAllocationMatrix() divides mix columns by that scale.
%   For a full-row-rank B, that is equivalent to B_norm = diag(scale) * B.

    mix_raw = pinv(B);
    scale = px4_allocation_scale_from_mix(mix_raw, normalize_rpy);
    D = diag(scale);
    B_norm = D * B;

    mix_norm = mix_raw;
    for axis = 1:min(numel(scale), size(mix_norm, 2))
        if scale(axis) > eps('single')
            mix_norm(:, axis) = mix_norm(:, axis) / scale(axis);
        end
    end

    mix_norm(abs(mix_norm) < 1e-3) = 0;
end

function scale = px4_allocation_scale_from_mix(mix, normalize_rpy)
%PX4_ALLOCATION_SCALE_FROM_MIX MATLAB version of PX4 updateControlAllocationMatrixScale().
%   mix is actuator_count x axis_count. Axis order follows PX4:
%   [roll pitch yaw thrust_x thrust_y thrust_z].

    axis_count = size(mix, 2);
    scale = ones(axis_count, 1);
    zero_tol = 1e-3;

    if normalize_rpy && axis_count >= 3
        num_non_zero_roll = nnz(abs(mix(:, 1)) > zero_tol);
        num_non_zero_pitch = nnz(abs(mix(:, 2)) > zero_tol);

        roll_norm_scale = 1;
        if num_non_zero_roll > 0
            roll_norm_scale = sqrt(sum(mix(:, 1).^2) / (num_non_zero_roll / 2));
        end

        pitch_norm_scale = 1;
        if num_non_zero_pitch > 0
            pitch_norm_scale = sqrt(sum(mix(:, 2).^2) / (num_non_zero_pitch / 2));
        end

        scale(1) = max(roll_norm_scale, pitch_norm_scale);
        scale(2) = scale(1);
        scale(3) = max(mix(:, 3));  % PX4 uses max(), not max(abs()).

    elseif axis_count >= 3
        scale(1:3) = 1;
    end

    if axis_count >= 6
        thrust_z = 6;
        scale(thrust_z) = 1;

        for axis = [6 5 4]
            col_abs = abs(mix(:, axis));
            num_non_zero_thrust = nnz(col_abs > eps('single'));

            if num_non_zero_thrust > 0
                scale(axis) = sum(col_abs) / num_non_zero_thrust;
            else
                scale(axis) = scale(thrust_z);
            end
        end
    end
end
