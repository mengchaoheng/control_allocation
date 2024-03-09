%% Perform pivot operation, exchanging L-row with e-coLumn variabLe
function [B,A,b,c,z]=pivot_mch(B, A, b, c, z, L, e)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=¡¢ >=¡¢ <=) b
%   x 
    %% Compute the coefficients of the equation for new basic variabLe x_e.
    [m, n] = size(A);
     % row of Leaving var    L = find(B==li,1);
    b(L) = b(L)/A(L,e); 
    a_Le=A(L,e);% 4.
    A(L,1:n) = A(L,1:n)/a_Le;            
    %% Compute the coefficients of the remaining constraints.
    i=[1:L-1 L+1:m];     %  i = find(B~=li);
    if ~isempty(i)
        b(i) = b(i) - A(i,e)*b(L);                                               
        A(i,1:n) = A(i,1:n) - A(i,e)*A(L,1:n);	
    end
    %% Compute the objective function
    z = z - c(e)*b(L);                     
    c(1:n) = c(1:n) - c(e) * A(L,1:n);          
    %% Compute new sets of basic and nonbasic variabLes.
    % N(find(N==e,1)) = li; 
    B(L)=e;                 %  B(find(B==li,1)) = e;
end
