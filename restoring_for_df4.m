l1=0.167;l2=0.069;k_v=3;
I_x=0.01149;
I_y=0.01153;
I_z=0.00487;
I=[I_x 0 0;0 I_y 0;0 0 I_z];
B=I\[-l1     0       l1     0;
     0      -l1     0       l1;
     l2    l2    l2    l2]*k_v;

[k,m] = size(B);
uMin=ones(m,1)*(-20)*pi/180;
uMax=ones(m,1)*20*pi/180;
mdes=[0;25;30];
global NumU
NumU=m;
LPmethod=2; % LPmethod should be an integer between 0 and 5. when LPmethod=2 set upper of lambda to Inf can't save this method!!! but big number is the same as that method based linprog
INDX=ones(1,m);  % active effectors
IN_MAT = [B     zeros(k,1)
          uMin' 0
          uMax' 0
          INDX  LPmethod];

% DP_LPCA vs DPscaled_LPCA
IN_MAT(1:3,end) = mdes;
u_DA1 = LPwrap(IN_MAT); % DP_LPCA




IN_MAT(end,end) =3; %DPscaled_LPCA 
u_DA2 = LPwrap(IN_MAT);


u_inv=min(max(pinv(B)*mdes, uMin), uMax)

% restoring
v_null=B*u_inv-B*u_inv;% =0

B_aug=[B;u_inv'];

v_aug= [zeros(k,1); -2];

u_null=pinv(B_aug)*v_aug

K_opt=2/(u_null'*u_null)

u_Pseudo = u_inv+K_opt*u_null % = pinv(B)*mdes

% update limits
uMin_new=uMin-u_inv;
uMax_new=uMax-u_inv;


% 0<K find K_max for K*u_null in limits
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
K_max
if K_max<K_opt
    u_DA_scaled = u_inv + K_max*u_null
else
    u_DA_scaled = u_inv + K_opt*u_null % u_Pseudo
end

u_DA_scaled'*u_DA_scaled-u_inv'*u_inv