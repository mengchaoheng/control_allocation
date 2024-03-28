#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include "allocator_dir_LPwrap_4.h"
#include "rt_nonfinite.h"
#include <set>
extern "C" {
    void allocator_dir_LPwrap_4(const float B[12], const float v[3],
                            const float umin[4], const float umax[4],
                            float u[4], float *z, unsigned int *iters);
}
Eigen::VectorXi setdiff(const Eigen::VectorXi& a, const Eigen::VectorXi& b) {
    // 将 Eigen::VectorXi 转换为 std::set<int>
    std::set<int> a_set(a.data(), a.data() + a.size());
    std::set<int> b_set(b.data(), b.data() + b.size());
    Eigen::VectorXi result;

    // 计算差集
    for (int elem : a_set) {
        if (b_set.find(elem) == b_set.end()) {
            result.conservativeResize(result.size() + 1);  // 动态调整结果向量的大小
            result(result.size() - 1) = elem;  // 将元素添加到结果向量中
        }
    }

    return result;
}
void BoundedRevisedSimplex(Eigen::MatrixXf A,
                            Eigen::RowVectorXf ct,
                            Eigen::VectorXf b,
                            Eigen::VectorXi& inB,
                            Eigen::VectorXf h,
                            Eigen::RowVector<bool, Eigen::Dynamic>& e,
                            int m, int n, int& itlim,
                            Eigen::VectorXf& y0,
                            int& errout) {
    // Solves the linear program:
    //      minimize c'y 
    //      subject to 
    //      Ay = b
    //      0<= y <= h
    float tol=1e-7;
    // std::cout << "Here is the float tol:\n" << tol << std::endl;

    // Eigen::VectorXi nind {{0, 1}}; // 1:n-m
    // 使用 LinSpaced() 创建等间隔序列，并转换为整数向量
    int n_m=n-m;
    Eigen::VectorXi nind = Eigen::VectorXi::LinSpaced(n_m, 0, n_m-1);
    std::cout << "Here is the VectorXi nind:\n" << nind << std::endl;

    // Eigen::VectorXi ind_all {{0, 1, 2, 3, 4}}; // 1:n
    Eigen::VectorXi ind_all = Eigen::VectorXi::LinSpaced(n, 0, n-1);
    std::cout << "Here is the VectorXi ind_all:\n" << ind_all << std::endl;
    
    Eigen::VectorXi inD = setdiff(ind_all, inB);
    std::cout << "Here is the VectorXi inD:\n" << inD << std::endl;

    // Initialize the matrix with replicated rows
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> e_matrix(m, n);
    e_matrix.rowwise() = e;
    // Output the initialized matrix
    std::cout << "Initialized e_matrix with every row the same:\n" << e_matrix << std::endl;
    A=e_matrix.select(A,-A);
    std::cout << "Here is the MatrixXf A:\n" << A << std::endl;
    ct=e.select(ct,-ct);
    std::cout << "Here is the RowVectorXf ct:\n" << ct << std::endl;
    
    Eigen::VectorXf h_reduce = e.transpose().select(0,h);
    std::cout << "Here is the VectorXf h_reduce:\n" << h_reduce << std::endl;
    Eigen::MatrixXf A_reduce = e_matrix.select(0,A);
    std::cout << "Here is the MatrixXf A_reduce:\n" << A_reduce << std::endl;

    Eigen::VectorXf Ah=A_reduce*h_reduce;
    std::cout << "Here is the VectorXf Ah:\n" << Ah << std::endl;
    std::cout << "Here is the before b:\n" << b << std::endl;
    b += Ah;
    std::cout << "Here is the after b:\n" << b << std::endl;
    // Initial Solution
    y0 = A(Eigen::all, inB).colPivHouseholderQr().solve(b); //https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    // std::cout << "A(Eigen::all, inB):\n" << A(Eigen::all, inB) << std::endl;
    std::cout << "The solution y0 is:\n" << y0 << std::endl;
    bool done = false;
    bool unbounded = false;
    
    // Main Simplex loop
    while ((!done || !unbounded) && (itlim > 0)) {
        itlim--;
        Eigen::VectorXf c=ct.transpose();
        Eigen::RowVectorXf lamt=A(Eigen::all, inB).transpose().colPivHouseholderQr().solve(c(inB)).transpose();
        // 
        std::cout << "Here is the RowVectorXf lamt:\n" << lamt << std::endl;
        Eigen::RowVectorXf rdt = ct(inD) - lamt*A(Eigen::all, inD);
        std::cout << "Here is the RowVectorXf rdt:\n" << rdt << std::endl;
        // Find minimum relative cost
        int qind;
        float minr = rdt.minCoeff(&qind);
        std::cout << "Here is the int qind:\n" << qind << std::endl;
        std::cout << "Here is the float minr:\n" << minr << std::endl;
        // If all relative costs are positive then the solution is optimal
        bool done = false;
        if (minr >= 0) {
            done = true;
            std::cout << "done"<< std::endl; //Check this condition
            break; // 退出循环
        }
        // Unknown to Enter the basis minimizes relative cost
        int qel=inD(qind);
        std::cout << "Here is the int qel:\n" << qel << std::endl;
        // Vector to enter in terms of the current Basis vector
        Eigen::VectorXf yq = A(Eigen::all, inB).colPivHouseholderQr().solve(A(Eigen::all, qel));
        std::cout << "Here is the VectorXf yq:\n" << yq << std::endl;
 
        if((yq.array().abs()<=tol).all())
        {
            unbounded = true;
            std::cout << "Solution is unbounded"<< std::endl; // Check this condition
            break;
        }
        // Compute ratio how much each current basic variable will have to move for the entering variable.
        Eigen::VectorXf rat = y0.array() / yq.array();
        std::cout << "Here is the VectorXf rat:\n" << rat << std::endl;
        // If yq < 0 then increasing variable when it leaves the basis will minimize cost
        Eigen::VectorXf hinB = h(inB);
        std::cout << "Here is the VectorXf hinB:\n" << hinB << std::endl;
        Eigen::Vector<bool, Eigen::Dynamic> indm = (yq.array() < 0).cast<bool>();
        rat.array() -=indm.select(hinB.array() / yq.array(),0);
        std::cout << "Here is the Vector indm:\n" << indm  << std::endl;
        std::cout << "Here is the VectorXf rat:\n" << rat  << std::endl;
        // If an element yq ~=0 then it doesn't change for the entering variable and shouldn't be chosen
        // one of the 3 version selcet: https://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html#a9e8e78c75887d4539071a0b7a61ca103
        Eigen::Vector<bool, Eigen::Dynamic> indz = (yq.array().abs() <= tol );
        rat= indz.select(std::numeric_limits<float>::infinity(), rat);
        std::cout << "Here is the Vector indz:\n" << indz  << std::endl;
        std::cout << "Here is the VectorXf rat:\n" << rat  << std::endl;
    
        // Variable to exit is moving to its minimum value
        int p;
        float minrat = rat.minCoeff(&p);
        std::cout << "Here is the int p:\n" << p  << std::endl;
        std::cout << "Here is the float minrat:\n" << minrat  << std::endl;

        // If the minimum ratio is zero, then the solution is degenerate and the entering
        // variable will not change the basis---invoke Bland's selection rule to avoid
        // cycling.
        if(std::abs(minrat)<=tol)
        {
            // Find negative relative cost
            std::cout << "Here is the rdt:\n" << rdt  << std::endl;
            Eigen::Vector<bool, Eigen::Dynamic> indm_nind = (rdt.array() < 0);
            std::cout << "Here is the Vector indm_nind:\n" << indm_nind  << std::endl;
            for (int i = 0; i < indm_nind.size(); ++i) {
                if (indm_nind(i)) {
                    qind = i;
                    std::cout << "Here is the int qind:\n" << qind  << std::endl;
                    break;
                }
            }
            
            
            qel=inD(qind);
            std::cout << "Here is the int qel:\n" << qel  << std::endl;
            yq = A(Eigen::all, inB).colPivHouseholderQr().solve(A(Eigen::all, qel));
            std::cout << "Here is the VectorXf yq:\n" << yq << std::endl;
            if((yq.array().abs()<=tol).all())
            {
                unbounded = true;
                std::cout << "Solution is unbounded"<< std::endl; //Check this condition
                break;
            }
            // Recompute rations and determine variable to leave
            rat = y0.array() / yq.array();
            std::cout << "Here is the VectorXf rat:\n" << rat  << std::endl;
            //If yq < 0 then increasing variable when it leaves the basis will minimize cost
            hinB = h(inB);
            std::cout << "Here is the VectorXf hinB:\n" << hinB << std::endl;
            indm = (yq.array() < 0).cast<bool>();
            std::cout << "Here is the Vector indm:\n" << indm  << std::endl;
            rat.array() -=indm.select(hinB.array() / yq.array(),0);
            std::cout << "Here is the VectorXf rat:\n" << rat  << std::endl;
            // If an element yq ~=0 then it doesn't change for the entering variable and shouldn't be chosen
            indz = (yq.array().abs() <= tol );
            std::cout << "Here is the Vector indz:\n" << indz  << std::endl;
            rat= indz.select(std::numeric_limits<float>::infinity(), rat);
            std::cout << "Here is the VectorXf rat:\n" << rat  << std::endl;
            // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
            minrat = rat.minCoeff(&p);
            std::cout << "Here is the int p:\n" << p  << std::endl;
            std::cout << "Here is the float minrat:\n" << minrat  << std::endl;
        }
        // Maintain the bounded simplex as only having lower bounds by recasting 
        // any variable that needs to move to its opposite bound.
        if (minrat >= h(qel)){ //Case 1: Entering variable goes to opposite bound and current basis is maintained
            e(qel) = !e(qel);
            A.col(qel) = -A.col(qel);
            b += A.col(qel) * h(qel);
            ct(qel) = -ct(qel);
            std::cout << " Case 1 "<< std::endl; 
            std::cout << "Here is the e:\n" << e  << std::endl;
            std::cout << "Here is the A:\n" << A  << std::endl;
            std::cout << "Here is the b:\n" << b  << std::endl;
            std::cout << "Here is the ct:\n" << ct  << std::endl;

        }else if(yq(p) > 0){ //Case 2: Leaving variable returns to lower bound (0)	
                int pel = inB(p);
                inB(p)= qel;
                inD(qind)= pel;
                std::cout << " Case 21 "<< std::endl; 
                std::cout << "Here is the pel:\n" << pel  << std::endl;
                std::cout << "Here is the inB:\n" << inB  << std::endl;
                std::cout << "Here is the inD:\n" << inD  << std::endl;
        }else{ //Case 2: Leaving variable moves to upper bound	
            std::cout << " Case 22 "<< std::endl; 
            int pel = inB(p);
            e(pel)=!e(pel);
            A.col(pel) = -A.col(pel);
            inB(p)= qel;
            inD(qind)= pel;
            ct(pel) = -ct(pel);
            std::cout << "Here is the b:\n" << b  << std::endl;
            b += A.col(pel) * h(pel);
            std::cout << "Here is the b:\n" << b  << std::endl;
            std::cout << "Here is the A.col(pel) * h(pel):\n" << A.col(pel) * h(pel)  << std::endl;
            
            std::cout << "Here is the pel:\n" << pel  << std::endl;
            std::cout << "Here is the e:\n" << e  << std::endl;
            std::cout << "Here is the A:\n" << A  << std::endl;
            std::cout << "Here is the inB:\n" << inB  << std::endl;
            std::cout << "Here is the inD:\n" << inD  << std::endl;
            std::cout << "Here is the ct:\n" << ct  << std::endl;
            
        }
        // Compute new Basic solution;
        y0 = A(Eigen::all, inB).colPivHouseholderQr().solve(b);
        std::cout << "Here is the y0:\n" << y0  << std::endl;
        // 在这里添加你的代码
        std::cout << " loop "<< std::endl; 
    }
    
    
    errout = unbounded;

}

void DP_LPCA(const Eigen::VectorXf& yd, const Eigen::MatrixXf& B, const Eigen::VectorXf& uMin, const Eigen::VectorXf& uMax, int& itlim,   
             const float& upper_lam, Eigen::VectorXf& u, int& errout)
{
    // DP_LPCA
    // Direction Preserving Control Allocation Linear Program

    // function [u, errout] = DP_LPCA(yd,B,uMin,uMax,itlim,upper_lam);

    // Solves the control allocation problem while preserving the
    // objective direction for unattainable commands. The solution
    // is found by solving the problem,
    // min -lambda,
    // s.t. B*u = lambda*yd, uMin<=u<=uMax, 0<= lambda <=1

    // For yd outside the AMS, the solution returned is that the
    // maximum in the direction of yd.

    // For yd strictly inside the AMS, the solution achieves
    // Bu=yd and m-n controls will be at their limits; but there
    // is no explicit preference to which solution will be 
    // returned. This limits the usefulness of this routine as
    // a practical allocator unless preferences for attainable solutions
    // are handled externally.

    // (For derivation of a similar formulation see A.1.2 and A.2.3 in the
    // text)


    // Inputs:
    //         yd [n]    = Desired objective
    //         B [n,m]   = Control Effectiveness matrix
    //         uMin[m,1] = Lower bound for controls
    //         uMax[m,1] = Upper bound for controls
    //         itlim     = Number of allowed iterations limit
    //                         (Sum of iterations in both branches)

    // Outputs:
    //         u[m,1]     = Control Solution
    //         errout     = Error Status code
    //                         0 = found solution
    //                         <0 = Error in finding initial basic feasible solution
    //                         >0 = Error in finding final solution
    //                         -1,1 = Solver error (unbounded solution)
    //                         -2   = Initial feasible solution not found
    //                         -3,3 = Iteration limit exceeded
    //         itlim      = Number of iterations remaining after solution found

    // Calls:
    //         BoundedRevisedSimplex = Bounded Revised Simplex solver (simplxuprevsol.m)

    // Notes:
    // If errout ~0 there was a problem in the solution. %

    // Error code < 0 implies an error in the initialization and there is no guarantee on
    // the quality of the output solution other than the control limits.
    // Error code > 0 for errors in final solution--B*u is in the correct direction and has
    // magnitude < yd, but B*u may not equal yd (for yd attainable)
    // or be maximized (for yd unattainable)

    // Modification History
    // 2024-03-28      MengChaoHeng  Original

    // % Reformula the direction-preserving control allocation
    // % problem to:
    // % min z=[0; -1]'[u; a]   s.t.  [B -v][u; a] = 0
    // %                                umin <= u <= umax
    // %                                   0 <= a
    // % and set x=u-umin, then
    // % min z=[0; -1]'[x; a]   s.t.  [B -v][x; a] = -B*umin
    // %                                0 <= x <= umax-umin
    // %                                0 <= a
    // B=[-0.4440         0    0.4440         0;
    //       0   -0.4440         0    0.4440;
    //    0.2070    0.2070    0.2070    0.2070]
    errout=0;
    int n=B.rows();
    int m=B.cols();

    if((yd.array().abs()<=std::numeric_limits<float>::epsilon()).all())
    {
        errout = -1;
        u.setZero();
        return;
    }

    Eigen::MatrixXf A(n,m+1);
    // yd*=-1;
    A << B, -yd;

    std::cout << "Here is the MatrixXf A:\n" << A << std::endl;
    Eigen::VectorXf b= -B * uMin;
    std::cout << "Here is the VectorXf b:\n" << b << std::endl;
    // Eigen::VectorXf c {{0.0, 0.0, 0.0, 0.0, -1.0}};
    Eigen::VectorXf zerosVector = Eigen::VectorXf::Zero(m);
    Eigen::VectorXf c(m + 1);
    c << zerosVector, -1.0f;
    std::cout << "Here is the VectorXf c:\n" << c << std::endl;
    Eigen::RowVectorXf ct=c.transpose();
    Eigen::VectorXf h(m+1);
    
    h  << uMax-uMin, upper_lam;
    std::cout << "Here is the VectorXf h:\n" << h << std::endl;

    // To find Feasible solution construct problem with appended slack variables
    // A.6.4 Initialization of the Simplex Algorithm of <Aircraft control allocation>
    // 计算 diag(sb)，即将 b 中大于 0 的元素设为 1，小于等于 0 的元素设为 -1
    Eigen::VectorXf sb = 2.0f * (b.array() > 0).cast<float>() - 1.0f;
    std::cout << "sb:\n" << sb << std::endl;
    std::cout << "sb.asDiagonal():\n" <<  Eigen::MatrixXf(sb.asDiagonal()) << std::endl;
    // 构建 Ai，先将 A 和 diag(sb) 水平拼接
    Eigen::MatrixXf Ai(A.rows(), A.cols() + n);
    Ai << A, Eigen::MatrixXf(sb.asDiagonal());
    std::cout << "Ai:\n" << Ai << std::endl;
    // 构建 ci
    Eigen::VectorXf ci(m + 1 + n);
    ci.head(m + 1).setZero(); // 设置前 m+1 个元素为零
    ci.tail(n).setOnes();     // 设置后 n 个元素为 1

    std::cout << "ci:\n" << ci << std::endl;

    // 构建 inBi
    Eigen::VectorXi inBi = Eigen::VectorXi::LinSpaced(n, m + 1, m + n);
    
    // 构建 ei
    Eigen::RowVector<bool, Eigen::Dynamic> ei(m + n + 1); 
    ei.setOnes(); // 设置所有元素为 1

    // 构建 hi
    Eigen::VectorXf hi(m + n + 1);
    hi.head(m+1) = h;
    hi.tail(n) = 2 * b.array().abs();

    // 输出 inBi、ei 和 hi
    std::cout << "inBi:\n" << inBi << std::endl;
    std::cout << "ei:\n" << ei << std::endl;
    std::cout << "hi:\n" << hi << std::endl;

    // Use Bounded Revised Simplex to find initial basic feasible point of
    // original program
    Eigen::VectorXf y1(n);
    y1.setZero();
    int errsimp=0;
    // [y1, inB1, e1,itlim, errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,m+n+1,itlim);
    BoundedRevisedSimplex(Ai, ci.transpose(), b, inBi, hi, ei, n, m+1+n, itlim, y1, errsimp);
    // heck that Feasible Solution was found
    if (itlim<=0)
    {   
        errout = -3;
        std::cout << " Too Many Iterations Finding Final Solution "<< std::endl; 
    }
    if((inBi.array()>m).any())
    {
        errout = -2;
        std::cout << " No Initial Feasible Solution found "<< std::endl; 
    }
    if(errsimp)
    {
        errout = -1;
        std::cout << " Solver error "<< std::endl; 
    }
    Eigen::VectorXf xout(m+1);
    xout.setZero();
    // Construct an incorrect solution to accompany error flags
    if(errout!=0)
    {
        // std::cout << "ToDo:\n"<< std::endl;
        Eigen::Vector<bool, Eigen::Dynamic> indv = (inBi.array() <= m).cast<bool>();
        for (int i = 0; i < indv.cols(); ++i) {
            if (indv(i)) {
                xout(inBi(i)) = y1(i);
            }
        }
        xout=ei(Eigen::seq(0,m)).transpose().select(xout,-xout+h);
    }
    else // No Error continue to solve problem
    {
        // Eigen::RowVector<bool, Eigen::Dynamic> e(m+1);
        // e << true, true, false, true, true;

        // Eigen::VectorXi inB(n); // inB 是一个整型向量
        // inB << 0, 1, 3;
        std::cout << "Here is the VectorXi inBi:\n" << inBi << std::endl;
        Eigen::VectorXf y0(n);
        y0.setZero();
        Eigen::RowVector<bool, Eigen::Dynamic> e2(m + 1); 
        e2=ei(Eigen::seq(0,m));
        std::cout << "Here is the RowVector e2:\n" << e2 << std::endl;
        BoundedRevisedSimplex(A, ct, b, inBi, h, e2, n, m+1, itlim, y0, errsimp);
        // outside
        std::cout << "itlim:\n" << itlim << std::endl;
        std::cout << "final e2:\n" << e2  << std::endl;
        std::cout << "final y0:\n" << y0  << std::endl;
        std::cout << "final inBi:\n" << inBi  << std::endl;
        
        xout(inBi)=y0;
        std::cout << "Here is the xout:\n" << xout  << std::endl;

        xout=e2.transpose().select(xout,-xout+h);
        std::cout << "reverse xout:\n" << xout  << std::endl;

        if (itlim<=0){
            errout = 3;
            std::cout << " Too Many Iterations Finding Final Solution "<< std::endl; 
        }

        if(errsimp){
            errout = 1;
            std::cout << " Solver error "<< std::endl; 
        }

    }

    u = xout(Eigen::seq(0,m-1))+uMin;

    if(xout(m)>1){ //Use upper_lam to prevent control surfaces from approaching position limits
        u /= xout(m);
    }
    std::cout << "finally u:\n" << u  << std::endl;
}

void DPscaled_LPCA(const Eigen::VectorXf& yd, const Eigen::MatrixXf& B, const Eigen::VectorXf& uMin, const Eigen::VectorXf& uMax, int& itlim, Eigen::VectorXf& u, int& errout)
{
    // This approach developed by Bodson (2001) is used to eliminate one of the constraints,
    // see A.5.3 Reduced Program Size without Secondary Optimization
    // Direction Preserving Control Allocation Linear Program
    //     Reduced formulation (Solution Scaled from Boundary)

    // function [u,itlim,errout] = DPscaled_LPCA(yd,B,uMin,uMax,itlim);

    // Solves the control allocation problem while preserving the
    // objective direction for unattainable commands. The reduced
    // dimension of the linear program passed to the Bounded Revised
    // Simplex solver is formed by forcing the solution to be on the
    // boundary of the AMS and eliminating the highest magnitude
    // objective by solving the other constraints in terms of it.

    // For yd outside the AMS, the solution returned is that the
    // maximum in the direction of yd
    // B*u= lamda*yd
    // max lamda s.t. uMin <= u <= uMax

    // Reducing the degrees of freedom elminates the problems of redundant
    // solutions for attainable objectives. If the desired objective is on the
    // interior of the AMS the solution is scaled from the solution on the
    // boundary, yielding the same controls as the Direct Allocation solution.

    // (In the text this solution is discussed in section A.5.3)

    // (See Bodson, M., "Evaluation of Optimization Methods for
    //         Control Allocation",  AIAA 2001-4223).

    // Inputs:
    //         yd [n]    = Desired objective
    //         B [n,m]   = Control Effectiveness matrix
    //         uMin[m,1] = Lower bound for controls
    //         uMax[m,1] = Upper bound for controls
    //         itlim     = Number of allowed iterations limit
    //                         (Sum of iterations in both branches)

    // Outputs:
    //         u[m,1]     = Control Solution
    //         errout     = Error Status code
    //                         0 = found solution
    //                         <0 = Error in finding initial basic feasible solution
    //                         >0 = Error in finding final solution
    //                         -1,1 = Solver error (unbounded solution)
    //                         -2   = Initial feasible solution not found
    //                         -3,3 = Iteration limit exceeded
    //         itlim      = Number of iterations remaining after solution found

    // Calls:
    //         BoundedRevisedSimplex = Bounded Revised Simplex solver (simplxuprevsol.m)

    // Notes:
    // If yd is close to zero, u = 0;

    // Error code < 0 implies an error in the initialization and there is no guarantee on
    // the quality of the output solution other than the control limits.
    // Error code > 0 for errors in final solution.

    // Modification History
    // 2024-03-28      MengChaoHeng  Original

    // Initialize error code to zero
    errout=0;
    int n=B.rows();
    int m=B.cols();
    
    // Locate the maximum magnitude element in the desired objective
    int iy;
    float my = yd.maxCoeff(&iy);
    if(my<=std::numeric_limits<float>::epsilon())
    {
        errout = -1;
        u.setZero();
        return;
    }
    Eigen::Vector<int, 1> ind_max_yd;
    ind_max_yd << iy; 
    std::cout << "Here is the ind_max_yd :\n" << ind_max_yd << std::endl;
    Eigen::VectorXi ind_all = Eigen::VectorXi::LinSpaced(n, 0, n-1);
    Eigen::VectorXi ind_other_yd = setdiff(ind_all, ind_max_yd);
    std::cout << "Here is the ind_other_yd:\n" << ind_other_yd << std::endl;
    Eigen::VectorXi ind_reorder(n);
    ind_reorder << ind_max_yd, ind_other_yd;
    std::cout << "Here is the ind_reorder:\n" << ind_reorder << std::endl;
    
    // 创建一个单位PermutationMatrix
    Eigen::PermutationMatrix<3, 3> pmat; // 3x3矩阵的PermutationMatrix
    pmat.setIdentity(); // 设置为单位矩阵
 
    // 将第iy列移动到第一列
    std::swap(pmat.indices()[0], pmat.indices()[iy]);
    std::cout << "Here is the B:\n" << B << std::endl;

    // 应用PermutationMatrix
    Eigen::MatrixXf Bt = pmat * B;
    std::cout << "Here is the Bt:\n" << Bt << std::endl;

    std::cout << "Here is the yd:\n" << yd << std::endl;
    Eigen::VectorXf ydt = pmat * yd;
    std::cout << "Here is the ydt:\n" << ydt << std::endl;

    // 交换第2行和第3行的内容
    std::swap(ydt(1), ydt(2));
    std::cout << "Here is the ydt:\n" << ydt << std::endl;
    Bt.row(1).swap(Bt.row(2));
    std::cout << "Here is the Bt:\n" << Bt << std::endl;

    // Convert into a LP problem
    Eigen::MatrixXf M(n-1,n);
    M.setZero();
    // Eigen::MatrixXf tmp=(-ydt(0)*Eigen::MatrixXf::Identity(n-1, n-1));

    M << ydt(Eigen::seq(1,n-1)), -ydt(0)*Eigen::MatrixXf::Identity(n-1, n-1);
    // std::cout << "Here is the tmp:\n" << tmp << std::endl;
    // std::cout << "Here is the ydt(Eigen::seq(1,n-1)):\n" << ydt(Eigen::seq(1,n-1)) << std::endl;

    std::cout << "Here is the M:\n" << M << std::endl;
    Eigen::MatrixXf A = M*Bt;
    Eigen::VectorXf b = -A*uMin;
    Eigen::VectorXf c = -Bt.transpose()*ydt;
    Eigen::VectorXf h = uMax-uMin;

    // Eigen::MatrixXf A(n,m+1);
    // yd*=-1;
    // A << B, -yd;

    std::cout << "Here is the MatrixXf A:\n" << A << std::endl;
    // Eigen::VectorXf b= -B * uMin;
    std::cout << "Here is the VectorXf b:\n" << b << std::endl;
    // Eigen::VectorXf c {{0.0, 0.0, 0.0, 0.0, -1.0}};
    // Eigen::VectorXf zerosVector = Eigen::VectorXf::Zero(m);
    // Eigen::VectorXf c(m + 1);
    // c << zerosVector, -1.0f;
    std::cout << "Here is the VectorXf c:\n" << c << std::endl;
    Eigen::RowVectorXf ct=c.transpose();
    // Eigen::VectorXf h(m+1);
    
    // h  << uMax-uMin, upper_lam;
    std::cout << "Here is the VectorXf h:\n" << h << std::endl;

    // To find Feasible solution construct problem with appended slack variables
    // A.6.4 Initialization of the Simplex Algorithm of <Aircraft control allocation>
    // 计算 diag(sb)，即将 b 中大于 0 的元素设为 1，小于等于 0 的元素设为 -1
    Eigen::VectorXf sb = 2.0f * (b.array() > 0).cast<float>() - 1.0f;
    std::cout << "sb:\n" << sb << std::endl;
    std::cout << "sb.asDiagonal():\n" <<  Eigen::MatrixXf(sb.asDiagonal()) << std::endl;
    // 构建 Ai，先将 A 和 diag(sb) 水平拼接
    Eigen::MatrixXf Ai(A.rows(), A.cols() + n);
    Ai << A, Eigen::MatrixXf(sb.asDiagonal());
    std::cout << "Ai:\n" << Ai << std::endl;
    // 构建 ci
    Eigen::VectorXf ci(m - 1 + n);
    ci.head(m).setZero(); // 设置前 m+1 个元素为零
    ci.tail(n-1).setOnes();     // 设置后 n 个元素为 1

    std::cout << "ci:\n" << ci << std::endl;

    // 构建 inBi
    Eigen::VectorXi inBi = Eigen::VectorXi::LinSpaced(n-1, m, m + n-2);
    
    // 构建 ei
    Eigen::RowVector<bool, Eigen::Dynamic> ei(m + n - 1); 
    ei.setOnes(); // 设置所有元素为 1

    // 构建 hi
    Eigen::VectorXf hi(m + n - 1);
    hi.head(m) = h;
    hi.tail(n-1) = 2 * b.array().abs();

    // 输出 inBi、ei 和 hi
    std::cout << "inBi:\n" << inBi << std::endl;
    std::cout << "ei:\n" << ei << std::endl;
    std::cout << "hi:\n" << hi << std::endl;

    // Use Bounded Revised Simplex to find initial basic feasible point of
    // original program
    Eigen::VectorXf y1(n-1);
    y1.setZero();
    int errsimp=0;
    // [y1, inB1, e1,itlim,errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n-1,m+n-1,itlim);
    BoundedRevisedSimplex(Ai, ci.transpose(), b, inBi, hi, ei, n-1, m-1+n, itlim, y1, errsimp);
    // heck that Feasible Solution was found
    if (itlim<=0)
    {   
        errout = -3;
        std::cout << " Too Many Iterations Finding Final Solution "<< std::endl; 
    }
    if((inBi.array()>m-1).any())
    {
        errout = -2;
        std::cout << " No Initial Feasible Solution found "<< std::endl; 
    }
    if(errsimp)
    {
        errout = -1;
        std::cout << " Solver error "<< std::endl; 
    }
    Eigen::VectorXf xout(m);
    xout.setZero();
    // Construct an incorrect solution to accompany error flags
    if(errout!=0)
    {
         // std::cout << "ToDo:\n"<< std::endl;
        xout(inBi(Eigen::seq(0,m-1))) = y1(Eigen::seq(0,m-1));
        xout=ei(Eigen::seq(0,m-1)).transpose().select(xout,-xout+h);

    }
    else // No Error continue to solve problem
    {
        Eigen::VectorXf y0(n-1);
        y0.setZero();
        Eigen::RowVector<bool, Eigen::Dynamic> e2(m); 
        e2=ei(Eigen::seq(0,m-1));
        BoundedRevisedSimplex(A, ct, b, inBi, h, e2, n-1, m, itlim, y0, errsimp);
        // outside
        std::cout << "itlim:\n" << itlim << std::endl;
        std::cout << "final e2:\n" << e2  << std::endl;
        std::cout << "final y0:\n" << y0  << std::endl;
        std::cout << "final inBi:\n" << inBi  << std::endl;
        // Construct solution to original LP problem from bounded simplex output
        // Set non-basic variables to 0 or h based on e2
        // Set basic variables to y2 or h-y2.
        
        xout(inBi)=y0;
        std::cout << "Here is the xout:\n" << xout  << std::endl;

        xout=e2.transpose().select(xout,-xout+h);
        std::cout << "reverse xout:\n" << xout  << std::endl;

        if (itlim<=0){
            errout = 3;
            std::cout << " Too Many Iterations Finding Final Solution "<< std::endl; 
        }

        if(errsimp){
            errout = 1;
            std::cout << " Solver error "<< std::endl; 
        }

    }
    // Transform Solution Back Into control variables
    // u(i) = x(i)+umin(i) if e(i)
    u = xout(Eigen::seq(0,m-1))+uMin;
    // Rescale controls so solution is not on boundary of Omega.
    Eigen::VectorXf rho = ydt.transpose()*Bt*u/(ydt.transpose()*ydt);
    std::cout << "rho:\n" << rho  << std::endl;
    if(rho(0)>1.0f){ 
        u /= rho(0);
    }
    std::cout << "finally u:\n" << u  << std::endl;
}


int main(int argc, char **argv)
{
    std::ifstream file("../../input.csv");
    std::vector<std::vector<double> > data;
    std::string line;
 
    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream sline(line);
            std::vector<double> row;
            std::string value;
 
            while (std::getline(sline, value, ',')) {
                row.push_back(std::stod(value));
            }
 
            data.push_back(row);
        }
        file.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
        return 1;
    }
 
    // 假设所有行的列数都相同，将第一行的列数作为数组的大小
    size_t columns = data[0].size();
    size_t num = data.size();
    double** array = new double*[data.size()];
    for (size_t i = 0; i < data.size(); ++i) {
        array[i] = new double[columns];
        for (size_t j = 0; j < columns; ++j) {
            array[i][j] = data[i][j];
        }
    }
    // 打开CSV文件进行写入
    std::ofstream outFile("../../output.csv");
    if (!outFile.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }
    
    // const float _B[3][4] = { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    // float B[12];
    // for (int i = 0; i < 3; i++)
    // {
    //     for(int j=0;j<4;j++)
    //     {
    //         B[i+3*j] = _B[i][j];
    //     }
    // }
    //    
    //      
    // 
    // float u_all[4]={ 0.0,  0.0,   0.0,   0.0};
    // size_t array_size = sizeof(u_all) / sizeof(u_all[0]);

    // for new implement using Eigen
    Eigen::MatrixXf B { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    int m=B.cols();

    std::cout << "Here is the matrix B:\n" << B << std::endl;
    Eigen::VectorXf uMin {{-0.3491, -0.3491, -0.3491, -0.3491}};
    Eigen::VectorXf uMax {{0.3491, 0.3491, 0.3491, 0.3491}};
    std::cout << "Here is the VectorXf uMin:\n" << uMin << std::endl;
    std::cout << "Here is the VectorXf uMax:\n" << uMax << std::endl;
    Eigen::VectorXf  yd(3);
    yd << -0.2064, 0.3253, 0.3187;
    std::cout << "Here is the MatrixXf yd:\n" << yd << std::endl;
    int itlim = 50; // 这里假设你的迭代次数初始值是50
    float upper_lam = 1e4f;
    int errout = 0;
    Eigen::VectorXf u(m);
    DP_LPCA(yd, B, uMin, uMax, itlim, upper_lam, u, errout);


    // for(int i=0;i<num;i++)
	// {
    //     float y_all[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};
    //     float z_all= 0.0;
    //     unsigned int iters_all= 0;
    //     float _uMin[4] ={};
    //     float _uMax[4] ={};
    //     for (int i = 0; i < 4; i++)
    //     {
    //         _uMin[i] =  -0.3491;
    //         _uMax[i] =  0.3491;
    //     }
    //     auto start = std::chrono::high_resolution_clock::now();
    //     allocator_dir_LPwrap_4(B, y_all, _uMin, _uMax, u_all, &z_all, &iters_all);
    //     auto finish = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> elapsed = finish - start;
    //     // std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
    //     // 使用循环打印数组元素
    //     // std::cout << "u_all: " << " ";
    //     // for (int i = 0; i < 4; ++i) {
    //     //     std::cout << u_all[i] << " ";
    //     // }
    //     // std::cout << std::endl; 

    //     // 写入CSV文件
    //     for (size_t i = 0; i < array_size; ++i) {
    //         outFile << u_all[i] << (i < array_size - 1 ? "," : "\n");
    //     }
    //     Eigen::MatrixXd mat(2, 2);
    // }

    // 关闭文件
    outFile.close();

    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;
    return 0;
}
