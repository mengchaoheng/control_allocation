
function [u_rest] = restoring(B,u,uMin,uMax)
% u have to be admissible
% restoring
tol=0.00001;
% by all(abs(null(B)'*u)) < eps or norm(null(B)'*u)<100*eps or rank([B_aug v_aug]) ~= rank(B_aug)
% for cpp is difficult to calc null(B) but we can calc
% norm(B*u_null)>0.00001 where u_null=pinv(B_aug)*v_aug;
if(norm(null(B)'*u)<tol) % a=0
    % null(B)
    % null(B)'
    % null(null(B)')
    % v_0=B*null(null(B)');
    % v_0
    % null(B)'*u
    u_rest=u;
    return;
end

[k,m] = size(B);
B_aug=[B;u'];

% update limits
uMin_new=uMin-u;
uMax_new=uMax-u;

a=-2;%a<0 => K>0  or  a>0 => K<0. so assume that K>0
v_aug= [zeros(k,1); a];

u_null=pinv(B_aug)*v_aug;

% R=rank(B_aug) = k
% B*u_null
% B_aug*u_null
K_opt=-a/(u_null'*u_null); % 

% u_Pseudo = u+K_opt*u_null % = pinv(B)*mdes

% 0<K find K_max for K*u_null in new limits
K_max=Inf;
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

if K_max<K_opt % or check if u + K_opt*u_null is admissible. If not, use K_max
    u_rest = u + K_max*u_null;
else
    u_rest = u + K_opt*u_null; % u_Pseudo
end

end