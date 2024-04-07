
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
    // Eigen::Vector<int, 1> ind_max_yd;
    // ind_max_yd << iy; 
    // std::cout << "Here is the ind_max_yd :\n" << ind_max_yd << std::endl;
    // Eigen::VectorXi ind_all = Eigen::VectorXi::LinSpaced(n, 0, n-1);
    // Eigen::VectorXi ind_other_yd = setdiff(ind_all, ind_max_yd);
    // std::cout << "Here is the ind_other_yd:\n" << ind_other_yd << std::endl;
    // Eigen::VectorXi ind_reorder(n);
    // ind_reorder << ind_max_yd, ind_other_yd;
    // std::cout << "Here is the ind_reorder:\n" << ind_reorder << std::endl;
    
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
    Eigen::MatrixXf Ai(A.rows(), A.cols() + n-1);
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
    float rho = ydt.dot(Bt*u)/(ydt.dot(ydt));
    std::cout << "rho:\n" << rho  << std::endl;
    if(rho>1.0f){ 
        u /= rho;
    }
    // std::cout << "finally u:\n" << u  << std::endl;
}
