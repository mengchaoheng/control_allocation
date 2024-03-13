%% Perform pivot operation, exchanging L-row with e-coLumn variabLe
function [B,A,b,C,z]=pivot_C(B, A, b, C, z, L, e,m,n)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=C*x   subj. to  A*x (=、 >=、 <=) b
%   x 
    %% Compute the coefficients of the equation for new basic variabLe x_e.
%     [m, n] = size(A);
     % row of Leaving var    L = find(B==li,1);
    b(L) = b(L)/A(L,e); 
    a_Le=A(L,e);% 4.
    A(L,1:n) = A(L,1:n)/a_Le;            
    %% Compute the coefficients of the remaining constraints.
    
%     i=[1:L-1 L+1:m];     %  i = find(B~=li);
    %     if ~isempty(i)
%         b(i) = b(i) - A(i,e)*b(L);                                               
%         A(i,1:n) = A(i,1:n) - A(i,e)*A(L,1:n);	
%     end
    b_L=b(L);
    for i=1:m
        if i~=L
            b(i) = b(i) - A(i,e)*b_L;
            a_je=A(i,e);
            for k=1:n
                A(i,k) = A(i,k) - a_je*A(L,k);
            end
        end
    end

    %% Compute the objective function
    C_e=C(e);
    z = z - C(e)*b(L);  
    for k=1:n
           C(k) = C(k) - C_e * A(L,k);
    end
%     C(1:n) = C(1:n) - C_e * A(L,1:n);      
    %% Compute new sets of basic and nonbasic variabLes.
    % N(find(N==e,1)) = li; 
    B(L)=e;                 %  B(find(B==li,1)) = e;
end
