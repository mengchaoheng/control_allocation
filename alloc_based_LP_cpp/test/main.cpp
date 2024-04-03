#include <iostream>
#include <matrix/math.hpp>
#include"ControlAllocation.h"
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
using namespace matrix;


int main() {
    std::ifstream file("/Users/mch/Proj/control_allocation/input.csv");
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
    std::ofstream outFile("/Users/mch/Proj/control_allocation/output.csv");
    if (!outFile.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }
    
    float _B[3][4] = { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
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
    float yd[3]={-0.430392439767736,-0.236610246030909,-0.0936906572928623};
    

    // 线性规划数据
    LinearProgrammingProblem<3, 5> problem;
    problem.tol=1e-7;
    problem.itlim = 10;
    //填数据
    problem.inB[0]=0;
    problem.inB[1]=1;
    problem.inB[2]=3;

    problem.e[0] = true;
    problem.e[1] = true;
    problem.e[2] = false;
    problem.e[3] = true;
    problem.e[4] = true;
    
    float upper_lam=1e4;
    for(int i=0; i<problem.m; ++i)
    {
        float temp=0;
        for(int j=0; j<problem.n-1; ++j)
        {
            problem.A[i][j] =_B[i][j];
            temp +=-_B[i][j]*_uMin[j];
        }
        problem.A[i][problem.n-1] =-yd[i];
        problem.b[i] = temp;
    }
    for(int i=0; i<problem.n-1; ++i)
    {
        problem.c[i] =0;
    }
    problem.c[problem.n-1] =-1;
    for(int i=0; i<problem.n; ++i)
    {
        problem.h[i] =_uMax[i]-_uMin[i];
    }
    problem.h[problem.n-1]=upper_lam;

    // 调用函数模板
    auto start = std::chrono::high_resolution_clock::now();
    LinearProgrammingResult<3, 5> result = BoundedRevisedSimplex(problem);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "BoundedRevisedSimplex execution time: " << elapsed.count() << "s\n";
    
    // 使用结果
    // result.y0, result.inB, result.e, result.errout
    int err = 0;
    float xout[problem.n];
    for(int i=0;i<problem.n;++i){
        xout[i]=0;
    }
    for(int i=0;i<problem.m;++i){
        xout[result.inB[i]]=result.y0[i];
    }
    for(int i=0;i<problem.n;++i){
        if(!result.e[i]){
            xout[i]=-xout[i]+problem.h[i];
        }
    }
    if(result.iters>=problem.itlim){
        err = 3;
        std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
    }
	if(result.errout)
    {
        err = 1;
        std::cout << "Solver error"<< std::endl;
    }
    float u_px4_matrix[problem.n-1];
    for(int i=0;i<problem.n-1;++i){
        u_px4_matrix[i]=xout[i]+_uMin[i];
    }
    if(xout[problem.n-1]>1){
        for(int i=0;i<problem.n-1;++i){
            u_px4_matrix[i]/=xout[problem.n-1];
        }
    }
    
    std::cout << "u_px4_matrix: [";
    for (size_t i = 0; i < problem.n-1; ++i) {
        std::cout << u_px4_matrix[i];
        if (i < problem.n-2) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;





    //飞机数据
    // 示例代码
    // 使用模板类时，可以指定 ControlSize 和 EffectSize 的具体值
    // 例如：
    float l1=0.148;float l2=0.069;float k_v=3;
    Aircraft<3, 4> df_4(_B, -0.3491, 0.3491); // 创建一个具有 4 个操纵向量和 3 个广义力矩的飞行器对象
    // 分配器数据：
    DP_LP_ControlAllocator<3, 4> DP_LPCA(df_4); // 创建一个控制分配器对象，用于具有 4 个操纵向量和 3 个广义力矩的飞行器(转化为线性规划问题，其维数和参数 <3, 4> 有关。)
    // 然后可以使用飞行器对象和控制分配器对象进行操作
   



    


    size_t array_size =4;
    // main loop
    for(int i=0;i<num;i++)
	{
        float yd[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};





        //=====================origin================================
        // bool flag=false;
        // for(int i=0;i<4;++i){
        //     if(std::abs(yd[i]) > problem.tol)
        //     {
        //         flag = true; // Check this condition
        //         break;
        //     }
        // }
        // if(!flag){
        //     for(int i=0;i<4;++i){
        //         u_px4_matrix[i]=0;
                
        //     }
        // }else{
        //     for(int i=0; i<problem.m; ++i)
        //     {
        //         problem.A[i][problem.n-1] =-yd[i];
        //     }
        //     start = std::chrono::high_resolution_clock::now();
        //     result = BoundedRevisedSimplex(problem);
        //     finish = std::chrono::high_resolution_clock::now();
        //     elapsed = finish - start;
        //     std::cout << "execution time: " << elapsed.count() << "s\n";
        //     // 使用结果
        //     // result.y0, result.inB, result.e, result.errout
        //     for(int i=0;i<problem.n;++i){
        //         xout[i]=0;
        //     }
        //     for(int i=0;i<problem.m;++i){
        //         xout[result.inB[i]]=result.y0[i];
        //     }
        //     for(int i=0;i<problem.n;++i){
        //         if(!result.e[i]){
        //             xout[i]=-xout[i]+problem.h[i];
        //         }
        //     }
        //     if(result.iters>=problem.itlim){
        //         err = 3;
        //         std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
        //     }
        //     if(result.errout)
        //     {
        //         err = 1;
        //         std::cout << "Solver error"<< std::endl;
        //     }
        //     for(int i=0;i<problem.n-1;++i){
        //         u_px4_matrix[i]=xout[i]+_uMin[i];
        //     }
        //     if(xout[problem.n-1]>1){
        //         for(int i=0;i<problem.n-1;++i){
        //             u_px4_matrix[i]/=xout[problem.n-1];
        //         }
        //     }

        // }
        // std::cout << "u_px4_matrix: [";
        // for (size_t i = 0; i < problem.n-1; ++i) {
        //     std::cout << u_px4_matrix[i];
        //     if (i < problem.n-2) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //=====================================================

        //=====================================================
         // 调用 u = DP_LP_ControlAllocator.allocateControl(yd)
        // float* u = new float[4];
        // start = std::chrono::high_resolution_clock::now();
        // u = DP_LPCA.allocateControl(yd);
        // // u = DP_LPCA.allocateControl(yd);
        // finish = std::chrono::high_resolution_clock::now();
        // elapsed = finish - start;
        // std::cout << "DP_LPCA.allocateControl execution time: " << elapsed.count() << "s\n";

        // std::cout << "u: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //=====================================================
        // float* u2 = new float[4];
        // start = std::chrono::high_resolution_clock::now();
        // u2 = DP_LPCA.allocateControl_bases_solver(yd);
        // finish = std::chrono::high_resolution_clock::now();
        // elapsed = finish - start;
        // std::cout << "DP_LPCA.allocateControl_bases_solver execution time: " << elapsed.count() << "s\n";

        // std::cout << "u2: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u2[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //=====================================================
        float* u3 = new float[4];
        start = std::chrono::high_resolution_clock::now();
        u3 = DP_LPCA.DP_LPCA(yd);
        // u = DP_LPCA.allocateControl(yd);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        std::cout << "DP_LPCA.DP_LPCA execution time: " << elapsed.count() << "s\n";
        std::cout << "u3: [";
        for (size_t i = 0; i < 4; ++i) {
            std::cout << u3[i];
            if (i < 3) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;


        // 写入CSV文件
        for (size_t i = 0; i < array_size; ++i) {
            outFile << u3[i] << (i < array_size - 1 ? "," : "\n");
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
