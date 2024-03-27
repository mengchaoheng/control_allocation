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

    
    const float _B[3][4] = { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    float B[12];
    for (int i = 0; i < 3; i++)
    {
        for(int j=0;j<4;j++)
        {
            B[i+3*j] = _B[i][j];
        }
    }
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

    float u_all[4]={ 0.0,  0.0,   0.0,   0.0};
    size_t array_size = sizeof(u_all) / sizeof(u_all[0]);
    // for new implement using Eigen
    Eigen::MatrixXf B_ac { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    // std::cout << "Here is the matrix B_ac:\n" << B_ac << std::endl;
    Eigen::VectorXf uMin_ac {{-0.3491, -0.3491, -0.3491, -0.3491}};
    Eigen::VectorXf uMax_ac {{0.3491, 0.3491, 0.3491, 0.3491}};
    // std::cout << "Here is the VectorXf uMin_ac:\n" << uMin_ac << std::endl;
    // std::cout << "Here is the VectorXf uMax_ac:\n" << uMax_ac << std::endl;
    Eigen::VectorXf  yd(3);
    yd << -0.2064, 0.3253, 0.3187;

    // std::cout << "Here is the MatrixXf yd:\n" << yd << std::endl;

    Eigen::MatrixXf A_ac(3,5);
    A_ac << B_ac, yd;
    std::cout << "Here is the MatrixXf A_ac:\n" << A_ac << std::endl;
    Eigen::VectorXf b_ac= -B_ac * uMin_ac;
    // std::cout << "Here is the VectorXf b_ac:\n" << b_ac << std::endl;
    Eigen::VectorXf c_ac {{0.0, 0.0, 0.0, 0.0, -1.0}};
    // std::cout << "Here is the VectorXf c_ac:\n" << c_ac << std::endl;
    Eigen::RowVectorXf ct_ac = c_ac.transpose();
    // std::cout << "Transposed vector of c_ac, RowVectorXf ct_ac:\n" << ct_ac << std::endl;


    Eigen::VectorXf h_ac(5);
    float upper_lam = 1e4f;
    h_ac  << uMax_ac-uMin_ac, upper_lam;
    // std::cout << "Here is the VectorXf h_ac:\n" << h_ac << std::endl;

    Eigen::VectorXi inB1_ac {{0, 1, 3}};
    // std::cout << "Here is the VectorXi inB1_ac:\n" << inB1_ac << std::endl;

    Eigen::VectorXi e1_ac(5);
    e1_ac << 1, 1, 0, 1, 1;
    // std::cout << "Here is the VectorXi e1_ac:\n" << e1_ac << std::endl;
    Eigen::RowVector<bool, Eigen::Dynamic> e_ac(5);
    e_ac << true, true, false, true, true;







    // // to BoundedRevisedSimplex
    int err=0;
    int m_A=3;
    int n_A=5;

    float tol=1e-10;
    // std::cout << "Here is the float tol:\n" << tol << std::endl;

    Eigen::VectorXi nind {{0, 1}}; // 1:n_A-m_A
    // std::cout << "Here is the VectorXi nind:\n" << nind << std::endl;

    Eigen::VectorXi ind_all {{0, 1, 2, 3, 4}}; // 1:n_A
    // std::cout << "Here is the VectorXi ind_all:\n" << ind_all << std::endl;
    Eigen::VectorXi inB(3); // inB 是一个整型向量
    inB << 0, 1, 3;
    // std::cout << "Here is the VectorXi inB:\n" << inB << std::endl;

    Eigen::VectorXi inD = setdiff(ind_all, inB);
    // std::cout << "Here is the VectorXi inD:\n" << inD << std::endl;

    // // Initialize the matrix with replicated rows
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> e_matrix(3, 5);
    e_matrix.rowwise() = e_ac;
    // Output the initialized matrix
    // std::cout << "Initialized e_matrix with every row the same:\n" << e_matrix << std::endl;
    A_ac=e_matrix.select(A_ac,-A_ac);
    // std::cout << "Here is the MatrixXf A_ac:\n" << A_ac << std::endl;

    ct_ac=e_ac.select(ct_ac,-ct_ac);
    std::cout << "Here is the RowVectorXf ct_ac:\n" << ct_ac << std::endl;
    
    Eigen::VectorXf h_ac_reduce = e_ac.transpose().select(0,h_ac);
    // std::cout << "Here is the VectorXf h_ac_reduce:\n" << h_ac_reduce << std::endl;
    Eigen::MatrixXf A_ac_reduce = e_matrix.select(0,A_ac);
    // std::cout << "Here is the MatrixXf A_ac_reduce:\n" << A_ac_reduce << std::endl;

    Eigen::VectorXf Ah=A_ac_reduce*h_ac_reduce;
    // std::cout << "Here is the VectorXf Ah:\n" << Ah << std::endl;
    b_ac += Ah;
    // std::cout << "Here is the VectorXf b_ac:\n" << b_ac << std::endl;
    // Initial Solution
    Eigen::VectorXf y0 = A_ac(Eigen::all, inB).colPivHouseholderQr().solve(b_ac); //https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    // std::cout << "A_ac(Eigen::all, inB):\n" << A_ac(Eigen::all, inB) << std::endl;
    // std::cout << "The solution y0 is:\n" << y0 << std::endl;
    bool done = false;
    bool unbounded = false;
    int itlim = 50; // 这里假设你的迭代次数初始值是100
    // Main Simplex loop
    while ((!done || !unbounded) && (itlim > 0)) {
        itlim--;
        Eigen::RowVectorXf lamt=A_ac(Eigen::all, inB).transpose().colPivHouseholderQr().solve(c_ac(inB)).transpose();
        // 



        std::cout << "Here is the RowVectorXf lamt:\n" << lamt << std::endl;
        Eigen::RowVectorXf rdt = ct_ac(inD) - lamt*A_ac(Eigen::all, inD);
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
        Eigen::VectorXf yq = A_ac(Eigen::all, inB).colPivHouseholderQr().solve(A_ac(Eigen::all, qel));
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
        Eigen::VectorXf hinB = h_ac(inB);
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
            yq = A_ac(Eigen::all, inB).colPivHouseholderQr().solve(A_ac(Eigen::all, qel));
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
        if (minrat >= h_ac(qel)){ //Case 1: Entering variable goes to opposite bound and current basis is maintained
            e_ac(qel) = !e_ac(qel);
            A_ac.col(qel) = -A_ac.col(qel);
            b_ac += A_ac.col(qel) * h_ac(qel);
            ct_ac(qel) = -ct_ac(qel);
            std::cout << " Case 1 "<< std::endl; 
            std::cout << "Here is the e_ac:\n" << e_ac  << std::endl;
            std::cout << "Here is the A_ac:\n" << A_ac  << std::endl;
            std::cout << "Here is the b_ac:\n" << b_ac  << std::endl;
            std::cout << "Here is the ct_ac:\n" << ct_ac  << std::endl;

        }else if(yq(p) > 0){ //Case 2: Leaving variable returns to lower bound (0)	
                int pel = inB(p);
                inB(p)= qel;
                inD(qind)= pel;
                std::cout << " Case 21 "<< std::endl; 
                std::cout << "Here is the pel:\n" << pel  << std::endl;
                std::cout << "Here is the inB:\n" << inB  << std::endl;
                std::cout << "Here is the inD:\n" << inD  << std::endl;
        }else{ //Case 2: Leaving variable moves to upper bound	
            int pel = inB(p);
            e_ac(pel)=!e_ac(pel);
            A_ac.col(pel) = -A_ac.col(pel);
            inB(p)= qel;
            inD(qind)= pel;
            ct_ac(pel) = -ct_ac(pel);
            std::cout << " Case 22 "<< std::endl; 
            std::cout << "Here is the pel:\n" << pel  << std::endl;
                
            std::cout << "Here is the e_ac:\n" << e_ac  << std::endl;
            std::cout << "Here is the A_ac:\n" << A_ac  << std::endl;
            std::cout << "Here is the inB:\n" << inB  << std::endl;
            std::cout << "Here is the inD:\n" << inD  << std::endl;
            std::cout << "Here is the ct_ac:\n" << ct_ac  << std::endl;
        }
        // Compute new Basic solution;
        y0 = A_ac(Eigen::all, inB).colPivHouseholderQr().solve(b_ac);
        std::cout << "Here is the y0:\n" << y0  << std::endl;
        // 在这里添加你的代码
        std::cout << " loop "<< std::endl; 
    }
    
    
    bool errout = unbounded;

    // outside
    std::cout << "itlim:\n" << itlim << std::endl;
    bool errsimp = errout;
    Eigen::RowVector<bool, Eigen::Dynamic> e_ac_out = e_ac;
    std::cout << "Here is the e_ac_out:\n" << e_ac_out  << std::endl;
    Eigen::VectorXf y_out =y0;
    std::cout << "Here is the y_out:\n" << y_out  << std::endl;
    Eigen::VectorXi inB_out=inB;
    std::cout << "Here is the inB_out:\n" << inB_out  << std::endl;
    Eigen::VectorXf xout(n_A);
    xout.setZero();
    xout(inB_out)=y_out;
    std::cout << "Here is the xout:\n" << xout  << std::endl;

    xout=e_ac_out.transpose().select(xout,-xout+h_ac);
    std::cout << "reverse xout:\n" << xout  << std::endl;


    if (itlim<=0){
        err = 3;
        std::cout << " Too Many Iterations Finding Final Solution "<< std::endl; 
    }

    if(errsimp){
        err = 1;
        std::cout << " Solver error "<< std::endl; 
    }
       
    Eigen::VectorXf u_simplxuprevsol = xout(Eigen::seq(0,n_A-2))+uMin_ac;
    std::cout << "reverse xout(Eigen::seq(0,n_A-2)):\n" << xout(Eigen::seq(0,n_A-2)) << std::endl;

    Eigen::VectorXf u_finally = u_simplxuprevsol;
    std::cout << "reverse u_simplxuprevsol:\n" << u_simplxuprevsol  << std::endl;
    if(xout(n_A-1)>1){ //Use upper_lam to prevent control surfaces from approaching position limits
        u_finally =u_simplxuprevsol /xout(n_A-1);
    }
    std::cout << "reverse u_finally:\n" << u_finally  << std::endl;









    






    



    for(int i=0;i<num;i++)
	{
        float y_all[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};
        float z_all= 0.0;
        unsigned int iters_all= 0;
        float _uMin[4] ={};
        float _uMax[4] ={};
        for (int i = 0; i < 4; i++)
        {
            _uMin[i] =  -0.3491;
            _uMax[i] =  0.3491;
        }
        auto start = std::chrono::high_resolution_clock::now();
        allocator_dir_LPwrap_4(B, y_all, _uMin, _uMax, u_all, &z_all, &iters_all);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        // std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
        // 使用循环打印数组元素
        // std::cout << "u_all: " << " ";
        // for (int i = 0; i < 4; ++i) {
        //     std::cout << u_all[i] << " ";
        // }
        // std::cout << std::endl; 

        // 写入CSV文件
        for (size_t i = 0; i < array_size; ++i) {
            outFile << u_all[i] << (i < array_size - 1 ? "," : "\n");
        }
        Eigen::MatrixXd mat(2, 2);


    }
    // 使用Eigen库
    // Eigen::MatrixXd mat(2, 2);
    // mat << 1, 2,
    //        3, 4;
    // std::cout << "Matrix mat:\n" << mat << std::endl;
    // 关闭文件
    outFile.close();

    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;
    return 0;
}
