#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include <chrono>
#include "allocator_dir_LPwrap_4.h"
#include "rt_nonfinite.h"
extern "C" {
    void allocator_dir_LPwrap_4(const float B[12], const float v[3],
                            const float umin[4], const float umax[4],
                            float u[4], float *z, unsigned int *iters);
}

#define AE_NO_EXCEPTIONS
using namespace alglib;

int main(int argc, char **argv)
{
    std::ifstream file("../../data.csv");
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
    real_2d_array A = "[[-0.4440, 0, 0.4440, 0, 0.2],[0, -0.4440, 0, 0.4440, 0.1], [0.2070, 0.2070, 0.2070, 0.2070, 0]]"; 
    real_1d_array AL = "[0,0,0]";
    real_1d_array AU = AL;
    real_1d_array C = "[0,0,0,0,-1]";
    real_1d_array S = "[0.05,0.05,0.05,0.05,2]";
    real_1d_array BndL = "[-0.3491,-0.3491,-0.3491,-0.3491,0]";
    real_1d_array BndU = "[+0.3491,+0.3491,+0.3491,+0.3491,+inf]";
    real_1d_array X;
    minlpstate STATE;
    minlpreport REP;
    double eps=1E-6;

    minlpcreate(5, STATE);
    minlpsetcost(STATE, C);
    minlpsetbc(STATE, BndL, BndU);
    minlpsetlc2dense(STATE, A, AL, AU, 3);
    minlpsetscale(STATE, S);
    minlpsetalgodss(STATE, eps);// faster

    // Solve
    auto start = std::chrono::high_resolution_clock::now();
    minlpoptimize(STATE);
    minlpresults(STATE, X, REP);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "firstly execution time: " << elapsed.count() << "s\n";

    float u_all[4]={ 0.0,  0.0,   0.0,   0.0};
    size_t array_size = sizeof(u_all) / sizeof(u_all[0]);


    for(int i=0;i<num;i++)
	{
        // Solve again
        // real_1d_array minus_v = "[0.2, 0.1, 0.1]"; 
        // -v
        A[0][4]=-data[i][0];
        A[1][4]=-data[i][1];
        A[2][4]=-data[i][2];
        minlpsetlc2dense(STATE, A, AL, AU, 3);

        start = std::chrono::high_resolution_clock::now();
        minlpoptimize(STATE);
        minlpresults(STATE, X, REP);
        finish = std::chrono::high_resolution_clock::now();
        // printf("X is: %s\n", X.tostring(4).c_str()); // 
        elapsed = finish - start;
        std::cout << "minlpoptimize execution time: " << elapsed.count() << "s\n";

        // vs matlab code
        double tmp[4] = {X(0), X(1), X(2), X(3)};
        if(X(4)>1) // ToDo: have to check the valid of the X 
        {
            for (int i = 0; i < 4; i++) 
            {
                tmp[i] /= X(4);
            }
        }
        alglib::real_1d_array delta;
        delta.setcontent(4,tmp);
        printf("delta is: %s\n", delta.tostring(4).c_str());
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
        start = std::chrono::high_resolution_clock::now();
        allocator_dir_LPwrap_4(B, y_all, _uMin, _uMax, u_all, &z_all, &iters_all);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        std::cout << "allocator_dir_LPwrap_4 execution time: " << elapsed.count() << "s\n";
        // 使用循环打印数组元素
        std::cout << "u_all: " << " ";
        for (int i = 0; i < 4; ++i) {
            std::cout << u_all[i] << " ";
        }
        std::cout << std::endl; 

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
