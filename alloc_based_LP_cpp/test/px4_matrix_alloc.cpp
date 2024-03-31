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

#include <matrix/math.hpp>

using namespace matrix;



/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */

extern "C" {
    void allocator_dir_LPwrap_4(const float B[12], const float v[3],
                            const float umin[4], const float umax[4],
                            float u[4], float *z, unsigned int *iters);
}

// 函数用于计算两个正整数集合的差
void setdiff(int setA[], int sizeA, int setB[], int sizeB, int result[], int& sizeResult) {
    sizeResult = 0;
    for (int i = 0; i < sizeA; ++i) {
        bool foundInB = false;
        // 检查当前 setA 中的元素是否在 setB 中
        for (int j = 0; j < sizeB; ++j) {
            if (setA[i] == setB[j]) {
                foundInB = true;
                break;
            }
        }
        // 如果当前元素不在 setB 中，则将其添加到结果中
        if (!foundInB) {
            result[sizeResult++] = setA[i];
        }
    }
}
int* generateSequence(int i, int n) {
    int* result = new int[n - i + 1]; // 动态分配数组内存

    for (int num = i, index = 0; num <= n; ++num, ++index) {
        result[index] = num;
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
    // //======================
    // // for new implement using Eigen
    // Eigen::MatrixXf B { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    // int m=B.cols();
    // std::cout << "Here is the matrix B:\n" << B << std::endl;
    // Eigen::VectorXf uMin {{-0.3491, -0.3491, -0.3491, -0.3491}};
    // Eigen::VectorXf uMax {{0.3491, 0.3491, 0.3491, 0.3491}};
    // std::cout << "Here is the VectorXf uMin:\n" << uMin << std::endl;
    // std::cout << "Here is the VectorXf uMax:\n" << uMax << std::endl;
    // Eigen::VectorXf  yd(3);
    // yd << 0.1, -0.1, 0.2;
    // std::cout << "Here is the MatrixXf yd:\n" << yd << std::endl;
    // int itlim = 500; // 这里假设你的迭代次数初始值是50
    // float upper_lam = 1e4f;
    // int errout = 0;
    // Eigen::VectorXf u(m);
    // u.setZero(); // important
    // // DP_LPCA(yd, B, uMin, uMax, itlim, upper_lam, u, errout);
    // DPscaled_LPCA(yd, B, uMin, uMax, itlim, u, errout);
    //===========================
    // for gen by matlab 
    const float _B[3][4] = { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    float _B_array[12];
    for (int i = 0; i < 3; i++)
    {
        for(int j=0;j<4;j++)
        {
            _B_array[i+3*j] = _B[i][j];
        }
    }
    float _uMin[4] ={};
    float _uMax[4] ={};
    for (int i = 0; i < 4; i++)
    {
        _uMin[i] =  -0.3491;
        _uMax[i] =  0.3491;
    }
    float u_all[4]={ 0.0,  0.0,   0.0,   0.0};
    size_t array_size = sizeof(u_all) / sizeof(u_all[0]);
    //===========================
    // for alloc_based_matrix
    // 注意 int& itlim，在多次调用时注意值的设置。
    // Solves the linear program:
    //      minimize c'y 
    //      subject to 
    //      Ay = b
    //      0<= y <= h
    auto start = std::chrono::high_resolution_clock::now();
    const int m=3;
    const int n=5;
    int inB[m]= {0, 1, 3};
    float tol=1e-7;
    int n_m=n-m;
    int* nind = generateSequence(0, n_m-1);
    int* ind_all = generateSequence(0, n-1);
    int inD[n_m]; 
    setdiff(ind_all, n, inB, m, inD, n_m);
    int itlim=100;
    float A[m][n];
    float b[m];
    float c[n];
    float h[n+1];
    float yd[m]{0.2, -0.1, 0.1};
    bool e[n]={true, true, false, true, true};;
    float upper_lam=1e4;
    for(int i=0; i<m; ++i)
    {
        float temp=0;
        for(int j=0; j<n-1; ++j)
        {
            A[i][j] =_B[i][j];
            temp +=-_B[i][j]*_uMin[j];
        }
        A[i][n-1] =-yd[i];
        b[i] = temp;
    }
    for(int i=0; i<n-1; ++i)
    {
        c[i] =0;
    }
    c[n-1] =-1;
    for(int i=0; i<n; ++i)
    {
        h[i] =_uMax[i]-_uMin[i];
    }
    h[n-1]=upper_lam;
    for(int i=0; i<m; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            if(!e[j])
            {
                A[i][j] *=-1;
                b[i]+=A[i][j]*h[j];
            }
        }
    }
    for(int j=0; j<n; ++j)
    {
        if(!e[j])
        {
            c[j] *=-1;
        }
    }
    //==============================
    matrix::SquareMatrix<float, 3> A_inB;
    matrix::Matrix<float, m, n-m> A_inD;
    matrix::Vector<float, m> c_inB;
    matrix::Vector<float, n-m> c_inD;
    for(int i=0; i<m; ++i)
    {
        for(int j=0; j<m; ++j)
        {
            A_inB(i,j)=A[i][inB[j]];
            if(j<n-m)
            {
                A_inD(i,j)=A[i][inD[j]];
            }
        }
        c_inB(i)=c[inB[i]];
    }
    for(int i=0; i<n-m; ++i)
    {
        c_inD(i)=c[inD[i]];
    }
    matrix::Vector<float, m> b_vec(b);
    //===============================
    //  %Initial Solution
    matrix::Vector<float, 3> y0 = inv(A_inB)*b_vec;
    bool done = false;
    bool unbounded = false;
    while ((!done  || !unbounded ) && (itlim > 0))
    {
        itlim = itlim-1;
        Matrix<float, 1UL, m> lamt= (inv(A_inB).transpose()*c_inB).transpose();
        Matrix<float, 1UL, n-m> rdt = c_inD.transpose()-lamt*A_inD;
        float minr;
        size_t qind;
        min(rdt.transpose(), minr, qind);
        if(minr >=0)  // If all relative costs are positive then the solution is optimal
        { 
            done = true;
            break;
        }
        int qel = inD[qind];  // Unknown to Enter the basis minimizes relative cost
        float vec_data[3] ={A[0][qel], A[1][qel], A[2][qel]};
        matrix::Vector<float, m> A_qel(vec_data);
        matrix::Vector<float, m> yq=inv(A_inB)* A_qel; // Vector to enter in terms of the current Basis vector
        bool flag=false;
        for(int i=0;i<m;++i){
            if(std::abs(yq(i)) > tol)
            {
                flag = true; // Check this condition
                break;
            }
        }
        if(!flag)
        {
            unbounded = true; // Check this condition
            break;
        }
        // Recompute rations and determine variable to leave
        matrix::Vector<float, m> rat;
        float hinB[m];
        for(int i=0;i<m;++i)
        {
            if(std::abs(yq(i))>tol)
            {
                rat(i)=y0(i)/yq(i);
                
            }
            else
            {
                rat(i)=INFINITY;
                /* code */
            }
        }
        for(int i=0;i<m;++i)
        {
            hinB[i]=h[inB[i]];
            if(yq(i)<0)
            {                    
                rat(i)-=hinB[i]/yq(i);
            }
        }
         // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        float minrat=rat(0);
        size_t p=0;
        min(rat, minrat, p);
        // If the minimum ratio is zero, then the solution is degenerate and the entering
        // variable will not change the basis---invoke Bland's selection rule to avoid
        // cycling.
        if (std::abs(minrat) <= tol)
        {
            //Find negative relative cost
            for(int i=0;i<n-m;++i)
            {
                if(rdt(1,i)<0){ //Note that since minr <0 indm is not empty   
                    qind=nind[i];
                    qel = inD[qind];//Unknown to Enter the basis is first indexed to avoid cycling
                    break;
                }
            }
            A_qel(0)=A[0][qel];
            A_qel(1)=A[1][qel];
            A_qel(2)=A[2][qel];
            yq=inv(A_inB)* A_qel;
            bool flag=false;
            for(int i=0;i<m;++i){
                if(std::abs(yq(i)) > tol)
                {
                    flag = true; // Check this condition
                    break;
                }
            }
            if(!flag)
            {
                unbounded = true; // Check this condition
                // break;
            }
            // Recompute rations and determine variable to leave
            // Recompute rations and determine variable to leave
            float hinB[m];
            for(int i=0;i<m;++i)
            {
                hinB[i]=h[inB[i]];
                if(std::abs(yq(i))>tol)
                {
                    rat(i)=y0(i)/yq(i);
                    if(yq(i)<0)
                    {                    
                        rat(i)-=hinB[i]/yq(i);
                    }
                }
                else
                {
                    rat(i)=INFINITY;
                    /* code */
                }
            }
            // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
            minrat=rat(0);
            p=0;
            min(rat, minrat, p);
        }
        if (minrat >= h[qel])
        {
            std::cout << " Case 1 "<< std::endl; 
            e[qel] =!e[qel];
            for(int i=0; i<m; ++i)
            {
                A[i][qel] *= -1;
                b_vec(i)+=A[i][qel]*h[qel];
            }
            c[qel] *= -1;

        }
        else if(yq(p) > 0)
        {
            std::cout << " Case 21 "<< std::endl; 
            int pel = inB[p];
            inB[p]= qel;
            inD[qind]= pel;
            // update x_inX
            for(int i=0; i<m; ++i)
            {
                A_inB(i,p)=A[i][qel];
            }
            for(int i=0; i<n-m; ++i)
            {
                c_inB(p)=c[qel];
            }
            for(int i=0; i<m; ++i)
            {
                A_inD(i,qind)=A[i][pel];
            }
            for(int i=0; i<n-m; ++i)
            {
                c_inD(qind)=c[pel];
            }
        }
        else
        {
            std::cout << " Case 22 "<< std::endl; 
            int pel = inB[p];
            e[pel]=!e[pel];
            for(int i=0; i<m; ++i)
            {
                A[i][pel] *= -1;
                b_vec(i)+=A[i][pel]*h[pel];
            }
            inB[p]= qel;
            inD[qind]= pel;
            c[pel] *= -1;
            // update x_inX
            for(int i=0; i<m; ++i)
            {
                A_inB(i,p)=A[i][qel];
            }
            
            for(int i=0; i<n-m; ++i)
            {
                c_inB(p)=c[qel];
            }
            for(int i=0; i<m; ++i)
            {
                A_inD(i,qind)=A[i][pel];
            }
            for(int i=0; i<n-m; ++i)
            {
                c_inD(qind)=c[pel];
            }
        }
        y0=inv(A_inB)* b_vec;
        std::cout << "y0:";
        y0.T().print();
    }
    bool errout = unbounded;  
    int err = 0;
    float xout[n];
    for(int i=0;i<n;++i){
        xout[i]=0;
    }
    for(int i=0;i<m;++i){
        xout[inB[i]]=y0(i);
    }
    for(int i=0;i<n;++i){
        if(!e[i]){
            xout[i]=-xout[i]+h[i];
        }
    }
    if(itlim<=0){
        err = 3;
        std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
    }
	if(errout)
    {
        err = 1;
        std::cout << "Solver error"<< std::endl;
    }
    float u_px4_matrix[n-1];
    for(int i=0;i<n-1;++i){
        u_px4_matrix[i]=xout[i]+_uMin[i];
    }
    if(xout[n-1]>1){
        for(int i=0;i<n-1;++i){
            u_px4_matrix[i]/=xout[n-1];
        }
    }

    std::cout << "u_px4_matrix: [";
    for (size_t i = 0; i < n-1; ++i) {
        std::cout << u_px4_matrix[i];
        if (i < n-2) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

	// main loop
    for(int i=0;i<1;i++)
	{
        float y_all[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};
        float z_all= 0.0;
        unsigned int iters_all= 0;
        // yd << (float) data[i][0],  (float) data[i][1],   (float) data[i][2];
        // std::cout << "yd:\n" << yd.transpose()  << std::endl;
        // auto start = std::chrono::high_resolution_clock::now();
        // allocator_dir_LPwrap_4(_B_array, y_all, _uMin, _uMax, u_all, &z_all, &iters_all);

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
        // 使用循环打印数组元素
        // std::cout << "u_all: " << " ";
        // for (int i = 0; i < 4; ++`i) {
        //     std::cout << u_all[i] << " ";
        // }
        // std::cout << std::endl; 
        // 写入CSV文件
        for (size_t i = 0; i < array_size; ++i) {
            outFile << u_all[i] << (i < array_size - 1 ? "," : "\n");
        }
    }
    // 关闭文件
    outFile.close();
    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;
    return 0;
}
