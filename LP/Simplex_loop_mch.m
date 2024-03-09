%% Iterate through simplex algorithm main loop
function [x,z]=Simplex_loop_mch(B, A, b, c, z)
% (c) mengchaoheng
% Last edited 2019-11
%   min z=c*x   subj. to  A*x (=�� >=�� <=) b
%   x 
    %% Initialization
    
    iters=0;
    [m,n] = size(A);
    while ~all(c>=0)                      % 3.~isempty(c(c(N)<0))
        e = find(c < 0, 1, 'first'); % ������������    % 4. e = N(find(c(N)<0,1))
        a_ie=A(:,e);
        ip=a_ie>0;
        delta=Inf*ones(m,1);
        if ~isempty(ip)
            delta(ip)=b(ip)./a_ie(ip);
        end
        [~,L]=min(delta);%ѡ����� ������
%         li = B(L);    % �����������               
        if delta(L) == Inf                         
            disp('System is unbounded!');
            return;  % 11.
        else
            [B,A,b,c,z] = pivot_mch(B,A,b,c,z,L,e); %�����������г����б任
        end
        iters = iters + 1;
    end
    x = zeros(n+m,1);                               % [13,16].
    x(B) = b; 
end
