#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "lp_tiny.h"
#include <chrono>

int main() {
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
 
    // 使用array...
    // 打印数据
    // for (const auto& row : data) {
    //     for (const auto& value : row) {
    //         std::cout << value << " ";
    //     }
    //     std::cout << std::endl;
    // }

    lp_tiny lp;
    lp_tiny_status status;
    double x[9];
    lp_tiny_init(&lp, 9, 7);
    
    lp.c[4] = -1;
    lp.b[2] = 0.2890;
    lp.b[3] = 0.6981;
    lp.b[4] = 0.6981;
    lp.b[5] = 0.6981;
    lp.b[6] = 0.6981;
    // i.e. A_{i,j} with n, m, is located at A[i+j*m]
    lp.A[0] = -0.4440;
    lp.A[2] = 0.2070;
    lp.A[3] = 1;

    lp.A[8] = -0.4440;
    lp.A[9] = 0.2070;
    lp.A[11] = 1;

    lp.A[14] = 0.4440;
    lp.A[16] = 0.2070;
    lp.A[19] = 1;

    lp.A[22] = 0.4440;
    lp.A[23] = 0.2070;
    lp.A[27] = 1;
    // -v
    lp.A[28] = 0;
    lp.A[29] = 0;
    lp.A[30] = 0;

    lp.A[38] = 1;
    lp.A[46] = 1;
    lp.A[54] = 1;
    lp.A[62] = 1;

    // 打开CSV文件进行写入
    std::ofstream outFile("../../output.csv");
    if (!outFile.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }
    double u[4];
    size_t array_size = sizeof(u) / sizeof(u[0]);

    printf("sphere test\n");
	for(int i=0;i<num;i++)
	{
        // -v
		lp.A[28] = -data[i][0];
		lp.A[29] = -data[i][1];
		lp.A[30] = -data[i][2];
        auto start = std::chrono::high_resolution_clock::now();
		lp_tiny_solve(&lp, x, NULL, NULL, &status, NULL);
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        std::cout << "firstly execution time: " << elapsed.count() << "s\n"; // 0.000196458s
		// printf(" Status: %s\n", lp_tiny_status_string(status));
		for(int i=0;i<4;i++)
		{
			u[i]=x[i]-0.3491;
		}
		if(x[4]>1)
		{
			for(int i=0;i<4;i++)
			{
				u[i]=u[i]/x[4];
			}
		}
		// printf(" Solution: %g, %g, %g, %g, f: %g\n", u[0], u[1], u[2], u[3], x[4]);
        // 写入CSV文件
        for (size_t i = 0; i < array_size; ++i) {
            outFile << u[i] << (i < array_size - 1 ? "," : "\n");
        }
        

	}
    // 关闭文件
    outFile.close();
    lp_tiny_destroy(&lp);
	return 0;



 
    // 释放内存
    for (size_t i = 0; i < data.size(); ++i) {
        delete[] array[i];
    }
    delete[] array;
    return 0;
}