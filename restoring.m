
function [u_rest] = restoring(B,u,uMin,uMax)
% u have to be admissible
% restoring
% by all(abs(null(B)'*u)) < eps or norm(null(B)'*u)<100*eps ?
if(norm(null(B)'*u)<100*eps) % a=0
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
K_max=1/eps;
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
    end
end

if K_max<K_opt
    u_rest = u + K_max*u_null;
else
    u_rest = u + K_opt*u_null; % u_Pseudo
end

end