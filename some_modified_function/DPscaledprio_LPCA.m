function [u,itlim,errout] = DPscaledprio_LPCA(yd,ye,B,uMin,uMax,itlim)
% Direction Preserving Control Allocation Linear Program
%     Reduced formulation (Solution Scaled from Boundary)
%
% function [u,itlim,errout] = DPscaled_LPCA(yd,B,uMin,uMax,itlim);
%
%    Solves the control allocation problem while preserving the
%  objective direction for unattainable commands. The reduced
%  dimension of the linear program passed to the Bounded Revised
%  Simplex solver is formed by forcing the solution to be on the
%  boundary of the AMS and eliminating the highest magnitude
%  objective by solving the other constraints in terms of it.
%
%  For yd outside the AMS, the solution returned is that the
%  maximum in the direction of yd
%    B*u= lamda*yd
%    max lamda s.t. uMin <= u <= uMax
%
%  Reducing the degrees of freedom elminates the problems of redundant
%  solutions for attainable objectives. If the desired objective is on the
%  interior of the AMS the solution is scaled from the solution on the
%  boundary, yielding the same controls as the Direct Allocation solution.
%  
%  (In the text this solution is discussed in section A.5.3)
%
%   (See Bodson, M., "Evaluation of Optimization Methods for
%          Control Allocation",  AIAA 2001-4223).
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    If yd is close to zero, u = 0;
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution.
%
% Modification History
%   2002      Roger Beck  Original ( DPcaLP2.m)
%   8/2014    Roger Beck  Update
[u_all,itlim_all,errout_all,rho] = DPscaled_LPCA(yd+ye,B,uMin,uMax,itlim);
if errout_all ~=0
    u=u_all;
    errout=errout_all;
    itlim=itlim_all;
else
    if rho>=1 %终端在AMS里，输出
        u=u_all;
        errout=errout_all;
        itlim=itlim_all;
    else %终端在AMS外
        [ue,itlime,erroute,rhoe] = DPscaled_LPCA(ye,B,uMin,uMax,itlim_all);
        if rhoe>1 %ye严格在AMS内，注意分配yd时标称位置已经变为ue
            uMin=uMin-ue;
            uMax=uMax-ue;
            [ud,itlimd,erroutd,rhod] = DPscaled_LPCA(yd,B,uMin,uMax,itlime); % rhod应大于1
            u=ud+ue; 
            errout=erroutd;
            itlim=itlimd;
        else  %ye严格在AMS上或在外,如果任意有限优先级，在这里递归，进一步去掉低优先级。
            u=ue;
            errout=erroute;
            itlim=itlime;
        end
    end
end
end
