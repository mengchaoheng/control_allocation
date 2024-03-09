function [u,error,errore,errord] = LPwrap(IN_MAT)
% Make single input, single output version of Linear Programming for use in
% Simulink via the MATLAB Fcn block
% IN_MAT = [B     yd         ye
%           umin' 0         Bu0(1)
%           umax' 0         Bu0(2)
%           INDX  LPmethod  Bu0(3)]
%
% 20140905  Created version to use Roger Beck's DB_LPCA program
% 20151206  Updated with Roger Beck's latest version of code and added 
%           LPmethod to select the various Linear Programming algorithms 

global NumU
% Get sizes
[k2,m1]=size(IN_MAT);
k=k2-3;
m=m1-2;
% If matrices too small, set contols to zero and return
if k<1 || m<1 || norm(IN_MAT)<1e-16
    u=zeros(NumU,1);
    return
end
% Partition input matrix into component matrices
B=IN_MAT(1:k,1:m);
v=IN_MAT(1:k,end-1);
ye=IN_MAT(1:k,end);
Bu0=IN_MAT(k+1:end,end);
umin=IN_MAT(k+1,1:m)';
umax=IN_MAT(k+2,1:m)';
LPmethod=IN_MAT(end,end-1);
% LPmethod should be an integer between 0 and 5
% 0 = DB_LPCA
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization Branch (1-norm)
        %      Control Error Minimization (1-norm)
% 1 = DBinf_LPCA
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization (1-norm)
        %      Control Error Minimization (inf-norm)
% 2 = DP_LPCA
        % Direction Preserving Control Allocation Linear Program    
% 3 =DPscaled_LPCA
        % Direction Preserving Control Allocation Linear Program
        %     Reduced formulation (Solution Scaled from Boundary)
% 4 = MO_LPCA
        % Mixed Optimization (Single Branch) Control Allocation Linear Program
        %    Objective Error Minimizing
        %    Control Error minimizing
% 5 = SB_LPCA
        % Single Branch Control Allocation Linear Program
        %    Direction Preserving
        %    Control Error minimizing
% If LPmethod is not one of the allowed options, set it to zero
if sum(LPmethod==[0 1 2 3 4 5 6 7 8 9 10 11 12 13])~=1
    LPmethod=0;
end

INDX=IN_MAT(end,1:m)';
m_act=sum(INDX);
% get active effectors
B_act=B(:,INDX>0.5); 
umin_act=umin(INDX>0.5);
umax_act=umax(INDX>0.5);

%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          wd [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution
%          wu [m,1]  = Weighting on control error
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%          lam       = Control Error Weighting (single parameter)
%          eMax[n,1] = Maximum Objective error
%          w [m,1]   = Control Error Weighting

% change variable names to be consistent with Roger's documentation
yd=v;
B=B_act;
uMin=umin_act;
uMax=umax_act;
% Some variables not specifically defined for generic case
% Values below are suggestions for now
% Users may add code here to specify alternative values
wd=ones(k,1);
up=zeros(m_act,1);
wu=0.1*ones(m_act,1);
emax=ones(k,1)*2e1;
itlim=100;
lam=0.1;
eMax=emax;
w=0.1*wu;

ue=zeros(4,1);
ud=zeros(4,1);

switch LPmethod
    case 0
        [u_act, feas, errout,itlim] = DB_LPCA(yd+ye-Bu0,B,wd,up,wu,emax,...
            uMin,uMax,itlim);
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization Branch (1-norm)
        %      Control Error Minimization (1-norm)
    case 1
        [u_act, feas, errout,itlim] = DBinf_LPCA(yd+ye-Bu0,B,wd,up,wu,emax,...
            uMin,uMax,itlim);
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization (1-norm)
        %      Control Error Minimization (inf-norm)
    case 2
        [u_act, errout] =  DP_LPCA(yd+ye-Bu0,B,uMin,uMax,itlim)
        % Direction Preserving Control Allocation Linear Program    
    case 3
        [u_act,itlim,errout,rho] = DPscaled_LPCA(yd+ye-Bu0,B,uMin,uMax,itlim);
        % Direction Preserving Control Allocation Linear Program
        %     Reduced formulation (Solution Scaled from Boundary)
    case 4
        [u_act,errout] = MO_LPCA(yd+ye-Bu0,B,up,lam, eMax,uMin,uMax,itlim);
        % Mixed Optimization (Single Branch) Control Allocation Linear Program
        %    Objective Error Minimizing
        %    Control Error minimizing
    case 5
        [u_act,errout] = SB_LPCA(yd+ye-Bu0,B,w,up,uMin,uMax,itlim);
        % Single Branch Control Allocation Linear Program
        %    Direction Preservingyd
        %    Control Error minimizing
    case 6
        [u_act,errout] = SBprio_LPCA(yd,ye,B,w,up,uMin,uMax,itlim);
        % Single Branch Control Allocation Linear Program
        %    Direction Preserving
        %    Control Error minimizing
    case 7
        [u_act,errout] = DPprio_LPCA(yd,ye,B,uMin,uMax,itlim);
        % Single Branch Control Allocation Linear Program
        %    Direction Preserving
        %    Control Error minimizing
    case 8
        [u_act,errout] = SBnew_LPCA(yd+ye-Bu0, B, w, up, uMin, uMax);
    case 9
        [u_act, errout] = DP3_LPCA(yd+ye-Bu0, B, uMin, uMax, itlim);
    case 10
        [u_act,errout,error,errore,errord] = DPnewprio_LPCA(yd, ye, B, uMin, uMax);
%         [u_act,itlim,errout] = DPscaledprio_LPCA(yd,ye,B,uMin,uMax,itlim);
    case 11
        [u_act,itlim,errout,rho] = DPscalednew_LPCA(yd, B, uMin, uMax, itlim);
    case 12
        [u_act,~,~] = dir_alloc_sim(yd, uMin, uMax, B);
    case 13
%         Ad=[B -yd -B yd; eye(7) -eye(7);-eye(7) eye(7)];
        Ad7=[-yd;0;0;0;0;0;0;1;0;0;0;0;0;0;-1];
        Ad14=[yd;0;0;0;0;0;0;-1;0;0;0;0;0;0;1];
        P_inv=[ -2    -4     6     0     0     0     0     0     0     0     0     0     0     0     0     0     0;
                 6     4    -6     0     0     0     0     0     0     0     0     0     0     0     0     0     0;
                -4     0     6     0     0     0     0     0     0     0     0     0     0     0     0     0     0;
                 2     4    -6     1     0     0     0     0     0     0     0     0     0     0     0     0     0;
                -6    -4     6     0     1     0     0     0     0     0     0     0     0     0     0     0     0;
                 4     0    -6     0     0     1     0     0     0     0     0     0     0     0     0     0     0;
                 0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0;
                 0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0;
                 0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0;
                 0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0;
                -2    -4     6     0     0     0     0     0     0     0     1     0     0     0     0     0     0;
                 6     4    -6     0     0     0     0     0     0     0     0     1     0     0     0     0     0;
                -4     0     6     0     0     0     0     0     0     0     0     0     1     0     0     0     0;
                 0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0;
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0;
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0;
                 0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1];
%                   A=[Ad_eye1 zeros(17,1) Ad_eye2 zeros(17,1) A1];
A =[ 1     0     0     1     2     2     0    -1     0     0    -1    -2    -2     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     1     0    -2    -3    -2     0     0    -1     0     2     3     2     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     1     2     2     1     0     0     0    -1    -2    -2    -1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0    -1    -2    -2     0     0     0     0     1     2     2     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     2     3     2     0     0     0     0    -2    -3    -2     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0    -2    -2    -1     0     0     0     0     2     2     1     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     1     0     0     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     1     0     0     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     1     0     0     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0;
     0     0     0     1     2     2     0     0     0     0    -1    -2    -2     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0;
     0     0     0    -2    -3    -2     0     0     0     0     2     3     2     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0;
     0     0     0     2     2     1     0     0     0     0    -2    -2    -1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0;
     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0;
     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0     0;
     0     0     0     0     0    -1     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     1];
        [u_act,~,~] = dir_alloc_six(uMin, uMax, A, Ad7,Ad14, P_inv);
end
u=zeros(NumU,1);
u(INDX>0.5,1)=u_act;
end