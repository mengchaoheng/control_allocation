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

#include "allocator_dir_LPwrap_4.h"
#include "rt_nonfinite.h"
extern "C" {
    void allocator_dir_LPwrap_4(const float B[12], const float v[3],
                            const float umin[4], const float umax[4],
                            float u[4], float *z, unsigned int *iters);
}


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
    

    //飞机数据
    // 示例代码
    // 使用模板类时，可以指定 ControlSize 和 EffectSize 的具体值
    // 例如：
    float l1=0.148;float l2=0.069;float k_v=3;
    Aircraft<3, 4> df_4(_B, -0.3491, 0.3491); // 创建一个具有 4 个操纵向量和 3 个广义力矩的飞行器对象
    DP_LP_ControlAllocator<3, 4> Allocator(df_4); // 创建一个控制分配器对象，用于具有 4 个操纵向量和 3 个广义力矩的飞行器(转化为线性规划问题，其维数和参数 <3, 4> 有关。)
    // 然后可以使用飞行器对象和控制分配器对象进行操作
   

    size_t array_size =4;

    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

    double total_elapsed1 = 0.0;
    double total_elapsed2 = 0.0;
    double total_elapsed3 = 0.0;
    double total_elapsed4 = 0.0;

    // main loop
    for(int i=0;i<num;i++)
	{
        float yd[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};
        // float yd[3]={0.1,  0.2,   -0.1}; // for test

        //==========================allocateControl===========================
        float u1[4]; int err1=0;
        start = std::chrono::high_resolution_clock::now();
        Allocator.allocateControl(yd, u1, err1); 
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed1 += elapsed.count();
        // std::cout << "Allocator.allocateControl execution time: " << elapsed.count() << "s\n";

        // std::cout << "u1: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u1[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //=========================DPscaled_LPCA=======have problem=====================INFO  [mixer_module] dir_alloc_sim time: 16
        float u2[4];int err2=0;float rho=0;
        start = std::chrono::high_resolution_clock::now();
        Allocator.DPscaled_LPCA(yd, u2, err2, rho);  
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed2 += elapsed.count();
        // std::cout << "DPscaled_LPCA rho: "<< rho <<std::endl; 
        // std::cout << "Allocator.DPscaled_LPCA execution time: " << elapsed.count() << "s\n";

        // std::cout << "u2: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u2[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //========================DP_LPCA=============================
        float u3[4];int err3=0;
        start = std::chrono::high_resolution_clock::now();
        Allocator.DP_LPCA(yd, u3, err3);  
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed3 += elapsed.count();
        // std::cout << "Allocator.DP_LPCA execution time: " << elapsed.count() << "s\n";

        // std::cout << "u3: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u3[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //========================allocator_dir_LPwrap_4 (generate by matlab) =============================
        float u4[4]={ 0.0,  0.0,   0.0,   0.0};
        float z_allocator_dir_LPwrap_4= 0.0;
        unsigned int iters_allocator_dir_LPwrap_4= 0;
        start = std::chrono::high_resolution_clock::now();
        allocator_dir_LPwrap_4(_B_array, yd, _uMin, _uMax, u4, &z_allocator_dir_LPwrap_4, &iters_allocator_dir_LPwrap_4); // allocator_dir_LPwrap_4 execution time: 7.08e-07s
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed4 += elapsed.count();
        // std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
        // // 使用循环打印数组元素
        // std::cout << "u4: " << " ";
        // for (int i = 0; i < 4; ++i) {
        //     std::cout << u4[i] << " ";
        // }
        // std::cout << std::endl; 

        // 写入CSV文件 change to u1 u2 u3 u4 for your test.
        for (size_t i = 0; i < array_size; ++i) {
            outFile << u4[i] << (i < array_size - 1 ? "," : "\n");
        }
    }
    // 求平均运行时间
    double average_elapsed1 = total_elapsed1 / num;
    std::cout << "Allocator.allocateControl Average execution time: " << average_elapsed1 << "s" << std::endl;
    double average_elapsed2 = total_elapsed2 / num;
    std::cout << "Allocator.DPscaled_LPCA Average execution time: " << average_elapsed2 << "s" << std::endl;
    double average_elapsed3 = total_elapsed3 / num;
    std::cout << "Allocator.DP_LPCA Average execution time: " << average_elapsed3 << "s" << std::endl;
    double average_elapsed4 = total_elapsed3 / num;
    std::cout << "allocator_dir_LPwrap_4 Average execution time: " << average_elapsed4 << "s" << std::endl;
    // running on M1 pro MacOS: 
    // Allocator.allocateControl Average execution time: 7.1138e-07s
    // Allocator.DPscaled_LPCA Average execution time: 5.03832e-07s
    // Allocator.DP_LPCA Average execution time: 1.50543e-06s
    // allocator_dir_LPwrap_4 Average execution time: 1.41764e-06s
    // 关闭文件
    outFile.close();
    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;








    

    return 0;
}
