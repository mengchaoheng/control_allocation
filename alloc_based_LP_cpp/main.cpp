#include <iostream>
#include <Eigen/Dense>

int main() {
    Eigen::MatrixXd mat(2, 2);
    mat << 1, 2,
           3, 4;
    
    std::cout << "Matrix mat:\n" << mat << std::endl;

    return 0;
}
