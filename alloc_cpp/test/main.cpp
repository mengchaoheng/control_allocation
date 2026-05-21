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
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

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
    // 打开CSV文件进行写入。所有 C++ 对比输出集中放在 results/cpp_outputs，
    // /results 已在 .gitignore 中，避免在工程根目录散落大 CSV 文件。
    const std::string outputRoot = "/Users/mch/Proj/control_allocation/results";
    const std::string outputDir = outputRoot + "/cpp_outputs";
    mkdir(outputRoot.c_str(), 0755);
    mkdir(outputDir.c_str(), 0755);
    const std::string outputPrefix = outputDir + "/";
    std::ofstream outFile(outputPrefix + "output.csv");
    std::ofstream outFile4(outputPrefix + "output_cpp_4.csv");
    std::ofstream outFile6(outputPrefix + "output_cpp_6.csv");
    std::ofstream outFile4DP(outputPrefix + "output_cpp_4_dp.csv");
    std::ofstream outFile4DPRaw(outputPrefix + "output_cpp_4_dp_raw.csv");
    std::ofstream outFile4DPscaled(outputPrefix + "output_cpp_4_dpscaled.csv");
    std::ofstream outFile4DPscaledRaw(outputPrefix + "output_cpp_4_dpscaled_raw.csv");
    std::ofstream outFile4Prio(outputPrefix + "output_cpp_4_prio.csv");
    std::ofstream outFile4PrioRaw(outputPrefix + "output_cpp_4_prio_raw.csv");
    std::ofstream outFile6DP(outputPrefix + "output_cpp_6_dp.csv");
    std::ofstream outFile6DPRaw(outputPrefix + "output_cpp_6_dp_raw.csv");
    std::ofstream outFile6DPscaled(outputPrefix + "output_cpp_6_dpscaled.csv");
    std::ofstream outFile6DPscaledRaw(outputPrefix + "output_cpp_6_dpscaled_raw.csv");
    std::ofstream outFile6Prio(outputPrefix + "output_cpp_6_prio.csv");
    std::ofstream outFile6PrioRaw(outputPrefix + "output_cpp_6_prio_raw.csv");
    if (!outFile.is_open() || !outFile4.is_open() || !outFile6.is_open()
        || !outFile4DP.is_open() || !outFile4DPRaw.is_open() || !outFile4DPscaled.is_open() || !outFile4DPscaledRaw.is_open()
        || !outFile4Prio.is_open() || !outFile4PrioRaw.is_open()
        || !outFile6DP.is_open() || !outFile6DPRaw.is_open() || !outFile6DPscaled.is_open() || !outFile6DPscaledRaw.is_open()
        || !outFile6Prio.is_open() || !outFile6PrioRaw.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }
    outFile << std::setprecision(17);
    outFile4 << std::setprecision(17);
    outFile6 << std::setprecision(17);
    outFile4DP << std::setprecision(17);
    outFile4DPRaw << std::setprecision(17);
    outFile4DPscaled << std::setprecision(17);
    outFile4DPscaledRaw << std::setprecision(17);
    outFile4Prio << std::setprecision(17);
    outFile4PrioRaw << std::setprecision(17);
    outFile6DP << std::setprecision(17);
    outFile6DPRaw << std::setprecision(17);
    outFile6DPscaled << std::setprecision(17);
    outFile6DPscaledRaw << std::setprecision(17);
    outFile6Prio << std::setprecision(17);
    outFile6PrioRaw << std::setprecision(17);
    
    float _I_x{0.01149f};//setting in the .sdf
	float _I_y{0.01153f};//setting in the .sdf
	float _I_z{0.00487f};//setting in the .sdf
	float _L_1{0.167f}; //setting in the .sdf
	float _L_2{0.069f}; //setting in the .sdf
	float _k{3.0f};	// USER_OMEGA_2_F, k  =_k_cv*_k_v*_k_v, setting k in the gazebo
    float _B[3][4] = {
        {-_L_1 / _I_x * _k, 0.0f, _L_1 / _I_x * _k, 0.0f},
        {0.0f, -_L_1 / _I_y * _k, 0.0f, _L_1 / _I_y * _k},
        {_L_2 / _I_z * _k, _L_2 / _I_z * _k, _L_2 / _I_z * _k, _L_2 / _I_z * _k}
    };


	matrix::Matrix<float, 4, 3> B_inv;
    B_inv.setZero();
	B_inv(0, 0)=-0.5f* _I_x/(_k*_L_1);
	B_inv(0, 2)=0.25f* _I_z/(_k*_L_2);

	B_inv(1, 1)=-0.5f* _I_y/(_k*_L_1);
	B_inv(1, 2)=0.25f* _I_z/(_k*_L_2);

	B_inv(2, 0)=0.5f* _I_x/(_k*_L_1);
	B_inv(2, 2)=0.25f* _I_z/(_k*_L_2);

	B_inv(3, 1)=0.5f* _I_y/(_k*_L_1);
	B_inv(3, 2)=0.25f* _I_z/(_k*_L_2);
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
    Aircraft<3, 4> df_4(_B, -0.3491f, 0.3491f); // 创建一个具有 4 个操纵向量和 3 个广义力矩的飞行器对象
    DP_LP_ControlAllocator<3, 4> Allocator(df_4); // 创建一个控制分配器对象，用于具有 4 个操纵向量和 3 个广义力矩的飞行器(转化为线性规划问题，其维数和参数 <3, 4> 有关。)
    // 然后可以使用飞行器对象和控制分配器对象进行操作
   

    // float _B_for6[3][6]    = { {-17.9880,   -8.9940,    8.9940,   17.9880,    8.9940,   -8.9940}, { 0,  -19.4726,  -19.4726,         0,   19.4726,   19.4726},{15.3231,   15.3231,   15.3231,   15.3231,   15.3231,   15.3231}};
    float l1_for6=0.2998f;float l2_for6=0.0664f;float k_for6=1.0f;
    float I_x_for6=0.05f; float I_y_for6=0.04f; float I_z_for6=0.013f;
    constexpr float DEG2RAD = 0.01745329251994329577f;
    const float d_for6 = 60.0f * DEG2RAD;
    float _B_for6[3][6] = {
        {-l1_for6/I_x_for6*k_for6, -l1_for6*std::cos(d_for6)/I_x_for6*k_for6,  l1_for6*std::cos(d_for6)/I_x_for6*k_for6,
          l1_for6/I_x_for6*k_for6,  l1_for6*std::cos(d_for6)/I_x_for6*k_for6, -l1_for6*std::cos(d_for6)/I_x_for6*k_for6},
        {0.0,  l1_for6*std::sin(d_for6)/I_y_for6*k_for6,  l1_for6*std::sin(d_for6)/I_y_for6*k_for6,
         0.0, -l1_for6*std::sin(d_for6)/I_y_for6*k_for6, -l1_for6*std::sin(d_for6)/I_y_for6*k_for6},
        {l2_for6/I_z_for6*k_for6, l2_for6/I_z_for6*k_for6, l2_for6/I_z_for6*k_for6,
         l2_for6/I_z_for6*k_for6, l2_for6/I_z_for6*k_for6, l2_for6/I_z_for6*k_for6}
    };
    constexpr float U_LIM = 40.0f * DEG2RAD;
    Aircraft<3, 6> df_6(_B_for6, -U_LIM, U_LIM);// 创建一个具有 6 个操纵向量和 3 个广义力矩的飞行器对象
    DP_LP_ControlAllocator<3, 6> Allocator_for6(df_6); // 创建一个控制分配器对象，用于具有 6 个操纵向量和 3 个广义力矩的飞行器(转化为线性规划问题，其维数和参数 <3, 6> 有关。)
    // 然后可以使用飞行器对象和控制分配器对象进行操作

    size_t array_size =4;

    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;

    double total_elapsed2 = 0.0;
    double total_elapsed3 = 0.0;
    double total_elapsed4 = 0.0;
    double total_elapsed5 = 0.0;
    double total_elapsed6 = 0.0;
    double total_elapsed7 = 0.0;
    double total_elapsed8 = 0.0;
    double total_elapsed9 = 0.0;
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

        //========================DP_LPCA=============================
        // MATLAB test.m counterpart:
        //   DP_LPCA(...), then restoring_cpp(...)
        // C++ allocator restoring() is aligned with restoring_cpp.m.
        float u2[4]; int err2=0; float rho2=0; float u2_raw[4];
        start = std::chrono::high_resolution_clock::now();
        Allocator.DP_LPCA(y_all, u2_raw, err2, rho2);
        Allocator.restoring(u2_raw, u2);
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
        float u3[4];int err3=0;float rho3=0;float u3_raw[4];
        start = std::chrono::high_resolution_clock::now();
        Allocator.DPscaled_LPCA(y_all, u3_raw, err3, rho3);
        Allocator.restoring(u3_raw, u3);
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
        // MATLAB test.m counterpart:
        //   DP_LPCA_prio(...), then restoring_cpp(...)
        // C++ DP_LPCA_prio() calls DP_LPCA_copy(), whose simplex
        // pivot/tie-break rules follow PCA/simplxuprevsol_tiebreak.m.
        float u4[4]; int err4=0; float rho4=0; float u4_raw[4];
        start = std::chrono::high_resolution_clock::now();
        Allocator.DP_LPCA_prio(m_higher, yd, u4_raw, err4, rho4);
        Allocator.restoring(u4_raw, u4);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed4 += elapsed.count();
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

        float m_higher_for6[3]={m_higher[0], m_higher[1], m_higher[2]};
        float yd_for6[3]={static_cast<float>(data[i][0]), static_cast<float>(data[i][1]), static_cast<float>(data[i][2])};
        float y_all_for6[3]={yd_for6[0]+m_higher_for6[0], yd_for6[1]+m_higher_for6[1], yd_for6[2]+m_higher_for6[2]};

        //========================DP_LPCA for 6 cs=============================
        // MATLAB test.m counterpart:
        //   DP_LPCA(...), then restoring_cpp(...)
        // C++ allocator restoring() is aligned with restoring_cpp.m.
        float u8[6]; int err8=0; float rho8=0; float u8_raw[6];
        start = std::chrono::high_resolution_clock::now();
        Allocator_for6.DP_LPCA(y_all_for6, u8_raw, err8, rho8);
        Allocator_for6.restoring(u8_raw, u8);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed8 += elapsed.count();

        //========================DPscaled_LPCA for 6 cs=============================
        float u9[6];int err9=0;float rho9=0;float u9_raw[6];
        start = std::chrono::high_resolution_clock::now();
        Allocator_for6.DPscaled_LPCA(y_all_for6, u9_raw, err9, rho9);
        Allocator_for6.restoring(u9_raw, u9);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed9 += elapsed.count();

        //========================DP_LPCA_prio for 6 cs=============================
        float u7[6]; int err7=0; float rho7=0; float u7_raw[6];
        start = std::chrono::high_resolution_clock::now();
        Allocator_for6.DP_LPCA_prio(m_higher_for6, yd_for6, u7_raw, err7, rho7);
        Allocator_for6.restoring(u7_raw, u7);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        total_elapsed7 += elapsed.count();

        //========================
        // 写入CSV文件：4舵和6舵分别输出，output.csv 保持为6舵 DP_LPCA 兼容输出。
        for (size_t i = 0; i < 4; ++i) {
            outFile4DP << u2[i] << (i < 4 - 1 ? "," : "\n");
            outFile4DPRaw << u2_raw[i] << (i < 4 - 1 ? "," : "\n");
            outFile4DPscaled << u3[i] << (i < 4 - 1 ? "," : "\n");
            outFile4DPscaledRaw << u3_raw[i] << (i < 4 - 1 ? "," : "\n");
            outFile4Prio << u4[i] << (i < 4 - 1 ? "," : "\n");
            outFile4PrioRaw << u4_raw[i] << (i < 4 - 1 ? "," : "\n");
        }
        for (size_t i = 0; i < 6; ++i) {
            outFile6DP << u8[i] << (i < 6 - 1 ? "," : "\n");
            outFile6DPRaw << u8_raw[i] << (i < 6 - 1 ? "," : "\n");
            outFile6DPscaled << u9[i] << (i < 6 - 1 ? "," : "\n");
            outFile6DPscaledRaw << u9_raw[i] << (i < 6 - 1 ? "," : "\n");
            outFile6Prio << u7[i] << (i < 6 - 1 ? "," : "\n");
            outFile6PrioRaw << u7_raw[i] << (i < 6 - 1 ? "," : "\n");
        }
        for (size_t i = 0; i < 4; ++i) {
            outFile4 << u2[i] << (i < 4 - 1 ? "," : "\n");
        }
        for (size_t i = 0; i < 6; ++i) {
            outFile6 << u8[i] << (i < 6 - 1 ? "," : "\n");
            outFile << u8[i] << (i < 6 - 1 ? "," : "\n");
        }
    }
    // 求平均运行时间
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
    double average_elapsed7 = total_elapsed7 / num;
    std::cout << "Allocator_for6.DP_LPCA_prio Average execution time: " << average_elapsed7 << "s" << std::endl;
    double average_elapsed8 = total_elapsed8 / num;
    std::cout << "Allocator_for6.DP_LPCA Average execution time: " << average_elapsed8 << "s" << std::endl;
    double average_elapsed9 = total_elapsed9 / num;
    std::cout << "Allocator_for6.DPscaled_LPCA Average execution time: " << average_elapsed9 << "s" << std::endl;
    // running on M1 pro MacOS: 
    // for float m_higher[3]={0.0,  0.0,  50.0}; //   have solution, m_higher attainable or m_higher+yd attainable
    // Allocator.DPscaled_LPCA Average execution time: 5.95756e-07s
    // Allocator.DP_LPCA Average execution time: 1.23543e-06s
    // allocator_dir_LPwrap_4 Average execution time: 1.23543e-06s
    // INV Average execution time: 1.82271e-08s
    // float m_higher[3]={0.0,  0.0,  150.0}; // No Initial Feasible Solution found m_higher+yd unattainable and m_higher unattainable
    // 关闭文件
    outFile.close();
    outFile4.close();
    outFile6.close();
    outFile4DP.close();
    outFile4DPRaw.close();
    outFile4DPscaled.close();
    outFile4DPscaledRaw.close();
    outFile4Prio.close();
    outFile4PrioRaw.close();
    outFile6DP.close();
    outFile6DPRaw.close();
    outFile6DPscaled.close();
    outFile6DPscaledRaw.close();
    outFile6Prio.close();
    outFile6PrioRaw.close();
    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;








    

    return 0;
}
