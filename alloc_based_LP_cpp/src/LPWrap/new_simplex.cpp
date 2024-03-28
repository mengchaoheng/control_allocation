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
void BoundedRevisedSimplex(Eigen::MatrixXf& A,
                            Eigen::RowVectorXf ct,
                            Eigen::VectorXf& b,
                            Eigen::VectorXi& inB,
                            Eigen::VectorXf& h,
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
    b += Ah;
    std::cout << "Here is the VectorXf b:\n" << b << std::endl;
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
    //         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)

    // Notes:
    // If errout ~0 there was a problem in the solution. %

    // Error code < 0 implies an error in the initialization and there is no guarantee on
    // the quality of the output solution other than the control limits.
    // Error code > 0 for errors in final solution--B*u is in the correct direction and has
    // magnitude < yd, but B*u may not equal yd (for yd attainable)
    // or be maximized (for yd unattainable)

    // Modification History
    // 2024-03-28      MengChaoHeng  Original
    errout=0;
    int n=B.rows();
    int m=B.cols();

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

    Eigen::RowVector<bool, Eigen::Dynamic> e(m+1);
    e << true, true, false, true, true;

    Eigen::VectorXi inB(n); // inB 是一个整型向量
    inB << 0, 1, 3;
    std::cout << "Here is the VectorXi inB:\n" << inB << std::endl;
    int errsimp=0;
    Eigen::VectorXf y0(n);
    y0.setZero();
    
    BoundedRevisedSimplex(A, ct, b, inB, h, e, n, m+1, itlim, y0, errsimp);
    // outside
    std::cout << "itlim:\n" << itlim << std::endl;
    std::cout << "final e:\n" << e  << std::endl;
    std::cout << "final y0:\n" << y0  << std::endl;
    std::cout << "final inB:\n" << inB  << std::endl;
    Eigen::VectorXf xout(m+1);
    xout.setZero();
    xout(inB)=y0;
    std::cout << "Here is the xout:\n" << xout  << std::endl;

    xout=e.transpose().select(xout,-xout+h);
    std::cout << "reverse xout:\n" << xout  << std::endl;

    if (itlim<=0){
        errout = 3;
        std::cout << " Too Many Iterations Finding Final Solution "<< std::endl; 
    }

    if(errsimp){
        errout = 1;
        std::cout << " Solver error "<< std::endl; 
    }
       
   u = xout(Eigen::seq(0,m-1))+uMin;

    if(xout(m)>1){ //Use upper_lam to prevent control surfaces from approaching position limits
        u /= xout(m);
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
    // %% but we use the Standard Forms for Linear Programming Problems
    // % min c'x subj. to A*x =b
    // %                  0 <= x
    // %% so we have to reformula the direction-preserving control allocation
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






    // 定义索引数量
    int nn = 3;

    // 使用 LinSpaced() 创建等间隔序列，并转换为整数向量
    Eigen::VectorXi index = Eigen::VectorXi::LinSpaced(nn, 1, nn);

    // 输出索引向量
    std::cout << "Indices vector:\n" << index.transpose() << std::endl;

    






    



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
