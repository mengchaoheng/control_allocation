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
    Eigen::MatrixXd mat(2, 2);
    mat << 1, 2,
           3, 4;
    std::cout << "Matrix mat:\n" << mat << std::endl;

    // 使用alglib库
    alglib::real_1d_array x = "[1,2,3]";
    alglib::real_1d_array y;
    
    return 0;
}
