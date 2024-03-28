#include <iostream>
#include <Eigen/Dense>
#include "stdafx.h"
#include <math.h>
#include "alglibinternal.h"
#include "alglibmisc.h"
#include "linalg.h"
#include "solvers.h"
#include "optimization.h"
#include "diffequations.h"
#include "specialfunctions.h"
#include "integration.h"
#include "statistics.h"
#include "interpolation.h"
#include "fasttransforms.h"
#include "dataanalysis.h"

int main() {
    // 使用Eigen库
    // Eigen::MatrixXd mat(2, 2);
    // mat << 1, 2,
    //        3, 4;
    // std::cout << "Matrix mat:\n" << mat << std::endl;

    // 使用alglib库
    alglib::real_1d_array x = "[1,2,3]";
    alglib::real_1d_array y;

//    Eigen::Matrix3f mat;
//     mat << 1, 2, 3,
//             4, 5, 6,
//             7, 8, 9;
 
//     // 创建一个单位PermutationMatrix
//     Eigen::PermutationMatrix<3, 3> pmat; // 3x3矩阵的PermutationMatrix
//     pmat.setIdentity(); // 设置为单位矩阵
 
//     // 将第二列移动到第一列
//     std::swap(pmat.indices()[0], pmat.indices()[1]);
 
//     // 应用PermutationMatrix
//     mat = pmat * mat;
 
//     std::cout << "Original matrix:\n" << mat << std::endl;
 
    int n = 5; // 假设n为矩阵的大小

    // 创建一个大小为n的单位矩阵，并乘以0.5
    Eigen::MatrixXd result = 0.5 * Eigen::MatrixXd::Identity(n, n);

    std::cout << "0.5 * eye(" << n << "):\n" << result << std::endl;


    
    return 0;
}
