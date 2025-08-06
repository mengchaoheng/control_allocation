
function [u_rest]  = restoring_opt(B,u,uMin,uMax)
% u have to be admissible, and  rank(B)==k
% restoring
% by all(abs(null(B)'*u)) < eps or norm(null(B)'*u)<100*eps or 
% rank([B_aug v_aug])=k+1 ~=rank(B_aug), or rank(B_aug)==k
% for cpp is difficult to calc null(B) but we can calc
% norm(B*pinv(B_aug)*v_aug)>0.00001 || all(abs(u) < eps) where  a=lager constant 
tol=0.00001;
[k,m] = size(B);
B_aug=[B;u'];
a=-2;
v_aug = [zeros(k,1); a];
u_null=pinv(B_aug)*v_aug;

% norm(B*u_null)>0.00001 <=> u'*u_null=0 && u!=0. 
% on the other hand u=0 => u_null=0 =>B*u_null=0
if( norm(B*u_null)>tol || all(abs(u) < tol)) 
% if (rank(B_aug)==k)
    % u_null
    % norm(B*u_null)
    K_opt=0;
    K_u=K_opt;
else
    K_opt=-a/(u_null'*u_null); % K_opt>0
    % 1.check whether u + K_opt*u_null is admissible. If not, use K_max
    % 2.or using the min one, but not rubostly (carfully for 0)
    % 3.or check K_opt*u_null, the same as aircraft control allocation
    % [isOuted, ~] = isOut(u + K_opt*u_null, uMin, uMax);
    % if isOuted  
        K_max=Inf;
        % update limits
        uMin_new=uMin-u;
        uMax_new=uMax-u;
        for i=1:m
            if(abs(u_null(i))<tol)
                continue;
            end
            if(u_null(i)>0)
                tmp=uMax_new(i)/u_null(i);
            else
                tmp=uMin_new(i)/u_null(i);
            end
            if(tmp<K_max) % find smaller
                K_max=tmp;
            end
        end
    if K_max < K_opt % or using the min one, but not rubostly
        K_u=K_max;
    else
        K_u=K_opt;
    end
end
u_rest = u + K_u*u_null;
end