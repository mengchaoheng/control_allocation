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
#include "wls_alloc_gen.h"
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
    
    // float _B[3][4] = { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
    float _B[3][4] = { {-43.6031,0.0,43.6031,0.0}, {0.0,-43.4519,0.0,43.4519},{42.5051,42.5051,42.5051,42.5051}};
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
    float I_array4[16];
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            I_array4[i + 4*j] = (i == j) ? 1.0f : 0.0f;
        }
    }
    float I_array3[9];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            I_array3[i + 3*j] = (i == j) ? 1.0f : 0.0f;
        }
    }


    //飞机数据
    // 示例代码
    // 使用模板类时，可以指定 ControlSize 和 EffectSize 的具体值
    // 例如：
    float l1=0.167;float l2=0.069;float k_v=3;
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
    double total_elapsed5 = 0.0;
    double total_elapsed6 = 0.0;
    // main loop
    for(int i=0;i<num;i++)
	{
        float m_higher[3]={0.0,  0.0,  0.}; // {0.0,  0.0,  40.};
        float yd[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};
        float y_all[3]={(float) data[i][0]+m_higher[0],  (float) data[i][1]+m_higher[1],   (float) data[i][2]+m_higher[2]};
        // float yd[3]={1.8729,  -3.2655,   0.1279}; // for test

        // change allocator parameters
        // if(i==1000) 
        // {
        //     for (int i = 0; i < 4; ++i) {
        //         Allocator.aircraft.upperLimits[i] = _uMax[i] + 0.2;
        //         Allocator.aircraft.lowerLimits[i] = _uMin[i] - 0.2;
        //          // update the limits
        //     }
        //     // update B
        //     for (int i = 0; i < 3; ++i) {
        //         for (int j = 0; j < 4; ++j) {
        //             Allocator.aircraft.controlEffectMatrix[i][j]= _B[i][j]*0.5;
        //         }
        //     }
        //     Allocator.isupdate = true; // have to set isupdate to true, otherwise the Allocator will not update the B and limits.
        // }else if(i==7000){  // need to reset the Allocator parameters
        //     for (int i = 0; i < 4; ++i) {
        //         Allocator.aircraft.upperLimits[i] = _uMax[i];
        //         Allocator.aircraft.lowerLimits[i] = _uMin[i];
        //          // update the limits
        //     }
        //     // update B
        //     for (int i = 0; i < 3; ++i) {
        //         for (int j = 0; j < 4; ++j) {
        //             Allocator.aircraft.controlEffectMatrix[i][j]= _B[i][j];
        //         }
        //     }
        //     Allocator.isupdate = true; // have to set isupdate to true, otherwise the Allocator will not update the B and limits.
        // }

        //==========================allocateControl===========================
        float u1[4]; int err1=0;
        start = std::chrono::high_resolution_clock::now();
        Allocator.allocateControl(y_all, u1, err1); 
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed1 += elapsed.count();
        // std::cout << "Allocator.allocateControl execution time: " << elapsed.count() << "s\n";
        // std::cout << "Allocator.err1: " << err1 << "\n";
        // std::cout << "u1: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u1[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //========================DP_LPCA=============================
        float u2[4];int err2=0;float rho2=0;float u2_tmp[4];
        start = std::chrono::high_resolution_clock::now();
        Allocator.DP_LPCA(y_all, u2_tmp, err2, rho2);  
        Allocator.restoring(u2_tmp,u2);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed2 += elapsed.count();
        // std::cout << "Allocator.DP_LPCA execution time: " << elapsed.count() << "s\n";
        // std::cout << "Allocator.err2: " << err2 << "\n";
        // std::cout << "u2: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u2[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //=========================DPscaled_LPCA=======have problem=====================INFO  [mixer_module] dir_alloc_sim time: 16
        float u3[4];int err3=0;float rho3=0;float u3_tmp[4];
        start = std::chrono::high_resolution_clock::now();
        Allocator.DPscaled_LPCA(y_all, u3_tmp, err3, rho3);  
        Allocator.restoring(u3_tmp,u3);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed3 += elapsed.count();
        // std::cout << "DPscaled_LPCA rho: "<< rho <<std::endl; 
        // std::cout << "Allocator.DPscaled_LPCA execution time: " << elapsed.count() << "s\n";
        // std::cout << "Allocator.err3: " << err3 << "\n";
        // std::cout << "u3: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u3[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //========================DP_LPCA_prio=============================
        float u4[4];int err4=0;float rho4=0;float u4_tmp[4]; 
	    // float m_lower[3]={30.0f,  0.0f,   -0.0f};
        int err41=0;float rho41=0;float u4_tmp1[4]; float m_tmp[3]={0.0,  0.0,  0.0f};
        start = std::chrono::high_resolution_clock::now();
        Allocator.DP_LPCA_copy(m_higher,yd, u4_tmp, err4, rho4); 
        if (err4<0){
            // std::cout << "Allocator.err4: " << err4 << "\n";
            Allocator.DP_LPCA_copy(m_tmp,m_higher, u4_tmp1, err41, rho41); 
            Allocator.restoring(u4_tmp1,u4);
        }else{
            Allocator.restoring(u4_tmp,u4);
            // std::cout << "Allocator.err4: " << err4 << "\n";
        }
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed4 += elapsed.count();
        // std::cout << "Allocator.DP_LPCA execution time: " << elapsed.count() << "s\n";
        // std::cout << "Allocator.err4: " << err4 << "\n";
        // std::cout << "Allocator.err41: " << err41 << "\n";
        // std::cout << "u4: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u4[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //========================allocator_dir_LPwrap_4 (generate by matlab) =============================
        float u5[4]={ 0.0,  0.0,   0.0,   0.0};
        float z_allocator_dir_LPwrap_4= 0.0;
        unsigned int iters_allocator_dir_LPwrap_4= 0;
        start = std::chrono::high_resolution_clock::now();
        allocator_dir_LPwrap_4(_B_array, y_all, _uMin, _uMax, u5, &z_allocator_dir_LPwrap_4, &iters_allocator_dir_LPwrap_4); // allocator_dir_LPwrap_4 execution time: 7.08e-07s
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed5 += elapsed.count();
        // std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
        // // 使用循环打印数组元素
        // std::cout << "u5: [";
        // for (size_t i = 0; i < 4; ++i) {
        //     std::cout << u5[i];
        //     if (i < 3) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        //=========================WLS_alloc_gen===========================
        float u6[4];float gam = 1e6f; float W0[4]={0.0f, 0.0f, 0.0f, 0.0f};  float u_d[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        start = std::chrono::high_resolution_clock::now();
        wls_alloc_gen(_B_array, y_all, _uMin, _uMax, I_array3, I_array4,
                   u_d, gam, u6, W0,
                   100, 4); 
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed6 += elapsed.count();
        // 写入CSV文件 change to u1 u2 u3 u4 for your test.
        for (size_t i = 0; i < array_size; ++i) {
            outFile << u4[i] << (i < array_size - 1 ? "," : "\n");
        }
    }
    // 求平均运行时间
    double average_elapsed1 = total_elapsed1 / num;
    std::cout << "Allocator.allocateControl Average execution time: " << average_elapsed1 << "s" << std::endl;
    double average_elapsed2 = total_elapsed2 / num;
    std::cout << "Allocator.DP_LPCA Average execution time: " << average_elapsed2 << "s" << std::endl;
    double average_elapsed3 = total_elapsed3 / num;
    std::cout << "Allocator.DPscaled_LPCA Average execution time: " << average_elapsed3 << "s" << std::endl;
    double average_elapsed4 = total_elapsed4 / num;
    std::cout << "Allocator.DP_LPCA_prio Average execution time: " << average_elapsed4 << "s" << std::endl;
    double average_elapsed5 = total_elapsed5 / num;
    std::cout << "allocator_dir_LPwrap_4 Average execution time: " << average_elapsed5 << "s" << std::endl;
    double average_elapsed6 = total_elapsed6 / num;
    std::cout << "wls_alloc_gen Average execution time: " << average_elapsed6 << "s" << std::endl;
    // running on M1 pro MacOS: 
    // for float m_higher[3]={0.0,  0.0,  50.0}; //   have solution, m_higher attainable or m_higher+yd attainable
    // Allocator.allocateControl Average execution time: 5.01194e-07s
    // Allocator.DPscaled_LPCA Average execution time: 5.95756e-07s
    // Allocator.DP_LPCA Average execution time: 1.23543e-06s
    // allocator_dir_LPwrap_4 Average execution time: 1.23543e-06s
    // Allocator.DP_LPCA_prio Average execution time: 1.66456e-06s
    // float m_higher[3]={0.0,  0.0,  150.0}; // No Initial Feasible Solution found m_higher+yd unattainable and m_higher unattainable
    // 关闭文件
    outFile.close();
    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;








    

    return 0;
}
