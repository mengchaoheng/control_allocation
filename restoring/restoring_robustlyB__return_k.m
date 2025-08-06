
function [u_rest,K_u,K_status,Deriv_ku_u] = restoring_robustlyB__return_k(B,u,uMin,uMax)
% u have to be admissible
% restoring
% by all(abs(null(B)'*u)) < eps or norm(null(B)'*u)<100*eps or rank([B_aug v_aug]) ~= rank(B_aug)
% for cpp is difficult to calc null(B) but we can calc 
% rank([B_aug v_aug]) ~= rank(B_aug)
[k,m] = size(B);
if(norm(null(B)'*u)<100*eps) % a=0
    % null(B)
    % null(B)'
    % null(null(B)')
    % v_0=B*null(null(B)');
    % v_0
    % null(B)'*u
    u_rest=u;
    % for test
    K_u=0;
    K_status=0; % k_u=k_opt=0
    Deriv_ku_u = zeros(m, m);
    return;
end


B_aug=[B;u'];
% rank(B_aug)
% update limits
uMin_new=uMin-u;
uMax_new=uMax-u;

a=-2;%a<0 => K>0  or  a>0 => K<0. so assume that K>0
v_aug= [zeros(k,1); a];

u_null=pinv(B_aug)*v_aug; %一定不为0

% R=rank(B_aug) = k
% B*u_null
% B_aug*u_null
K_opt=-a/(u_null'*u_null); % 

% u_Pseudo = u+K_opt*u_null % = pinv(B)*mdes

%
Deriv_uNull_ud = zeros(m, m);
P=B_aug*B_aug';
for i=1:m
    E=zeros(k+1,m);
    E(end,i)=1;
    Deriv_uNull_ud(:,i)=(-B_aug'*pinv(P)*(B_aug*E'+E*B_aug')*pinv(P)+E'*pinv(P))*v_aug;
end


% 0<K find K_max for K*u_null in new limits
K_max=Inf;
for i=1:m

    if(u_null(i)>0)
        if(abs(u_null(i))<eps)
            u_null(i)=eps;
        end
        tmp=uMax_new(i)/u_null(i);
    else
        if(abs(u_null(i))<eps)
            u_null(i)=-eps;
        end
        tmp=uMin_new(i)/u_null(i);
    end
    if(tmp<K_max) % find smaller
        K_max=tmp;
        max_index=i;
        if(u_null(i)>0)
            Deriv_k_u=( -u_null(i) - Deriv_uNull_ud(i,i)*uMax_new(i) ) / u_null(i)^2;
        else
            Deriv_k_u=( -u_null(i) - Deriv_uNull_ud(i,i)*uMin_new(i) ) / u_null(i)^2;
        end
    end
end
% 每次都运行u_rest=min(max(u + K_opt*u_null, uMin), uMax);会导致没有充分利用可达集

% [is_saturated, saturated_mask] = isSaturated(u + K_opt*u_null, uMin, uMax);
[isOuted, ~] = isOut(u + K_opt*u_null, uMin, uMax);
if ~isOuted % ((K_opt <=  K_max) || ~isOuted )     % 到这里，k_opt一定不为 0
    u_rest = u + K_opt*u_null; % u_Pseudo
    % for test
    K_u=K_opt;
    if K_opt==0
        K_max 
    end
    K_status=1;
    % for Deriv in inv
    Deriv_k_uNull=2*a*u_null'/(u_null'*u_null)^2;

    Deriv_ku_u=u_null*Deriv_k_uNull*Deriv_uNull_ud+K_opt*Deriv_uNull_ud;
    
    
else %K_max<K_opt  % k_max可能为0. 
    u_rest = u + K_max*u_null;
    % for test
    
    K_u=K_max;
    K_status=2;
    % for Deriv in AS but not in inv
    Deriv_k_ud=zeros(1,m);
    % if K_max==0
    %     Deriv_ku_u=zeros(m, m);
    % else
    %     Deriv_k_ud(max_index)=Deriv_k_u;
    %     Deriv_ku_u= u_null*Deriv_k_ud +K_max*Deriv_uNull_ud;
    % end
    Deriv_k_ud(max_index)=Deriv_k_u;
    Deriv_ku_u= u_null*Deriv_k_ud +K_max*Deriv_uNull_ud;
    if K_max~=0
        % u_null*Deriv_k_ud
        % u_null
    end

end

