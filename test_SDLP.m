clc;
clear all;
close all;
% If we use the Standard Forms for Linear Programming Problems
% min c'x subj. to A*x <= b
%                    0 <= x
%% so we have to reformula the direction-preserving control allocation
% problem to:
% min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
%                                umin <= u <= umax
%                                   0 <= a
% and set x=u-umin, then
% min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
%                                        x <= umax-umin
%                                        0 <= x 
%                                        0 <= a
% min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] <= -B*umin
%                              [-B v][x; a] <=  B*umin
%                              [I  0][x; a] <= umax-umin
%                                        -x <= 0 
%                                        -a <= 0
% set X=[x;a]
l1=0.148;l2=0.069;k_v=3;
B=k_v*[-l1     0       l1     0;
     0      -l1     0       l1;
     l2    l2    l2    l2];
[k,m] = size(B);
umin=ones(m,1)*(-20)*pi/180;
umax=ones(m,1)*20*pi/180;
v=[0.2;0.1;0.1];
A=[B -v;-B v; eye(m) zeros(m,1)];
b=[-B*umin; B*umin; umax-umin];
c=[zeros(m,1); -1];
[X,fval,exitflag,output,lambda]= linprog(c,A,b,[],[],zeros(1+m,1),[]);
u = X(1:m)+umin;
a = X(m+1);
% Scale down u if a>1
if a>1
    u = u/a;
end
u

%% for SDLP 
% A b c is 
% A=[B -v;-B v; eye(m) zeros(m,1);-eye(m) zeros(m,1);zeros(1,m) -1]
% b=[-B*umin; B*umin; umax-umin;0;0;0;0;0]
% c=[zeros(m,1); -1]
% or
A
b
c

% setup A b and c to sdlp_example and then run :
% ➜  build git:(main) ✗ make          
% [ 50%] Building CXX object CMakeFiles/sdlp_example.dir/example/sdlp_example.cpp.o
% [100%] Linking CXX executable sdlp_example
% [100%] Built target sdlp_example
% ➜  build git:(main) ✗ ./sdlp_example
% optimal sol: 0.155327 0.426713   0.6981   0.6981  1.20496
% u: -0.193773 0.0776133     0.349     0.349
% optimal obj: -1.20496
u_SDLP=[-0.193773; 0.0776133;     0.349;     0.349];
% optimal obj is the same, but B*u_SDLP != v, u != u_SDLP   
B*u_SDLP

