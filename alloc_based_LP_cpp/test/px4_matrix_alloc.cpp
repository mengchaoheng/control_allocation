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

int main(int argc, char **argv)
{

    // test data
    float data1[] = {1,2,3,4,5};
    float data2[] = {6,7,8,9,10};
    Vector<float, 5> v1(data1);
    Vector<float, 5> v2(data2);

    v2.print();
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
    for(int i=0;i<num;i++)
	{
        float y_all[3]={(float) data[i][0],  (float) data[i][1],   (float) data[i][2]};
        float z_all= 0.0;
        unsigned int iters_all= 0;
        // yd << (float) data[i][0],  (float) data[i][1],   (float) data[i][2];
        // std::cout << "yd:\n" << yd.transpose()  << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        allocator_dir_LPwrap_4(_B_array, y_all, _uMin, _uMax, u_all, &z_all, &iters_all);

        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        // std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
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
