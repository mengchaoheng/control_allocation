#include <iostream>

#include "sdlp/sdlp.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    // int m = 2 * 7;
    int m = 10;
    Eigen::Matrix<double, 5, 1> x;        // decision variables
    Eigen::Matrix<double, 5, 1> c;        // objective coefficients
    Eigen::Matrix<double, -1, 5> A(m, 5); // constraint matrix
    Eigen::VectorXd b(m);                 // constraint bound
    A <<   -0.4440,         0,              0.4440,    0,        -0.2000,
            0,             -0.4440,         0,         0.4440,   -0.1000,
            0.2070,         0.2070,         0.2070,    0.2070,   -0.1000,
            0.4440,         0,             -0.4440,    0,         0.2000,
            0,              0.4440,         0,        -0.4440,    0.1000,
           -0.2070,        -0.2070,        -0.2070,   -0.2070,    0.1000,
            1.0000,         0,              0,         0,         0,
            0,              1.0000,         0,         0,         0,
            0,              0,              1.0000,    0,         0,
            0,              0,              0,         1.0000,    0;


    b <<0, 0, 0.2890, 0, 0, -0.2890, 0.6981, 0.6981, 0.6981, 0.6981;


    c << 0.0, 0.0, 0.0, 0.0, -1.0;



    // c << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    // A << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //     0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    //     -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //     0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0,
    //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0;
    // b << 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0;

    double minobj = sdlp::linprog<5>(c, A, b, x);

    // std::cout << "prob:\n"
    //           << std::endl;
    // std::cout << "     min x1 + ... + x7," << std::endl;
    // std::cout << "     s.t. x1 <=  6,  x2 <=  5, ..., x7 <= 0," << std::endl;
    // std::cout << "          x1 >= -1,  x2 >= -2,  ..., x7 >= -7.\n"
    //           << std::endl;

    Eigen::Matrix<double, 4, 1> u_={x(0)-0.3491, x(1)-0.3491, x(2)-0.3491, x(3)-0.3491};
    Eigen::Matrix<double, 4, 1> u=u_;
    if(minobj>=1)
    {
        for(int i=0;i<4;i++)
        {
            u(i)=u_(i)/minobj;
        }
    }
    std::cout << "optimal sol: " << x.transpose() << std::endl;
    std::cout << "u: " << u.transpose() << std::endl;
    std::cout << "optimal obj: " << minobj << std::endl;

    return 0;
}