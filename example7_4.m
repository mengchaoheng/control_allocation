B=[0.7073 -0.7073 -3.4956 -3.0013 3.0013 3.4956 2.1103  
    1.1204 1.1204 -0.7919 -1.2614 -1.2614 -0.7919 0.0035  
    -0.3309 0.3309 -0.1507 -0.3088 0.3088 0.1507 -1.2680 ];
uMin = [  -0.9599  -0.9599  -0.5236  -0.5236  -0.5236  -0.5236  -0.5236]';
uMax =  [ 0.4363  0.4363  0.5236  0.5236  0.5236  0.5236  0.5236  ]';
mdes =[ -1.8832  1.4055  1.0133 ]';
itlim=5e2;
[k,m] = size(B);
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
u_DA1 = LPwrap(IN_MAT) % DP_LPCA




IN_MAT(end,end) =3; %DPscaled_LPCA 
u_DA2 = LPwrap(IN_MAT);


v_null=B*u_DA1-B*u_DA2; % =0

B_aug=[B;u_DA2'];

v_aug= [zeros(k,1); -2];

u_null=pinv(B_aug)*v_aug;

K_opt=2/(u_null'*u_null);

u_Pseudo = u_DA2+K_opt*u_null; % = pinv(B)*mdes

% update limits
uMin_new=uMin-u_DA2;
uMax_new=uMax-u_DA2;

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
if K_max<K_opt
    u_DA_scaled = u_DA2 + K_max*u_null
else
    u_DA_scaled = u_DA2 + K_opt*u_null % u_Pseudo
end

% Comparing u_Pseudo with uMin and uMax
% 0<K find K_max for K*u_Pseudo in limits
K=1/eps;
for i=1:m

    if(u_Pseudo(i)>0)
        if(abs(u_Pseudo(i))<eps)
            u_Pseudo(i)=eps;
        end
        tmp=uMax(i)/u_Pseudo(i);
    else
        if(abs(u_Pseudo(i))<eps)
            u_Pseudo(i)=-eps;
        end
        tmp=uMin(i)/u_Pseudo(i);
    end
    if(tmp<K) % find smaller
        K=tmp;
    end
end
u_Pseudo_scaled=u_Pseudo* K
B*u_Pseudo_scaled

% B*u_Pseudo_scaled is defferent from mdes=B*u_DA_scaled
B*u_DA_scaled

B*u_DA_scaled-B*u_DA2



u_DA_scaled'*u_DA_scaled-u_DA2'*u_DA2