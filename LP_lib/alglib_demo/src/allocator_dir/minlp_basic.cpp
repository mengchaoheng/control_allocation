#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include <chrono>
#include "dir_alloc_sim.h"
#include <iostream>
#define AE_NO_EXCEPTIONS
using namespace alglib;

int main(int argc, char **argv)
{
    // try
    // {
    //
    // This example demonstrates how to minimize
    //
    //     F(x0,x1) = -0.1*x0 - x1
    //
    // subject to box constraints
    //
    //     -1 <= x0,x1 <= +1 
    //
    // and general linear constraints
    //
    //     x0 - x1 >= -1
    //     x0 + x1 <=  1
    //
    // We use dual simplex solver provided by ALGLIB for this task. Box
    // constraints are specified by means of constraint vectors bndl and
    // bndu (we have bndl<=x<=bndu). General linear constraints are
    // specified as AL<=A*x<=AU, with AL/AU being 2x1 vectors and A being
    // 2x2 matrix.
    //
    // NOTE: some/all components of AL/AU can be +-INF, same applies to
    //       bndl/bndu. You can also have AL[I]=AU[i] (as well as
    //       BndL[i]=BndU[i]).
    //
        // real_2d_array a = "[[1,-1],[1,+1]]";
        // real_1d_array al = "[-1,-inf]";
        // real_1d_array au = "[+inf,+1]";
        // real_1d_array c = "[-0.1,-1]";
        // real_1d_array s = "[1,1]";
        // real_1d_array bndl = "[-1,-1]";
        // real_1d_array bndu = "[+1,+1]";
        // real_1d_array x;
        // minlpstate state;
        // minlpreport rep;

        // minlpcreate(2, state);

    //
    // Set cost vector, box constraints, general linear constraints.
    //
    // Box constraints can be set in one call to minlpsetbc() or minlpsetbcall()
    // (latter sets same constraints for all variables and accepts two scalars
    // instead of two vectors).
    //
    // General linear constraints can be specified in several ways:
    // * minlpsetlc2dense() - accepts dense 2D array as input; sometimes this
    //   approach is more convenient, although less memory-efficient.
    // * minlpsetlc2() - accepts sparse matrix as input
    // * minlpaddlc2dense() - appends one row to the current set of constraints;
    //   row being appended is specified as dense vector
    // * minlpaddlc2() - appends one row to the current set of constraints;
    //   row being appended is specified as sparse set of elements
    // Independently from specific function being used, LP solver uses sparse
    // storage format for internal representation of constraints.
    //
        // minlpsetcost(state, c);
        // minlpsetbc(state, bndl, bndu);
        // minlpsetlc2dense(state, a, al, au, 2);

    //
    // Set scale of the parameters.
    //
    // It is strongly recommended that you set scale of your variables.
    // Knowing their scales is essential for evaluation of stopping criteria
    // and for preconditioning of the algorithm steps.
    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
    //
        // minlpsetscale(state, s);

    // Solve
        // minlpoptimize(state);
        // minlpresults(state, x, rep);
        // printf("%s\n", x.tostring(3).c_str()); // EXPECTED: [0,1]


    // % (c) mengchaoheng
    // % Last edited 2024-03
    // %   min z=c*x   subj. to  A*x (=、 >=、 <=) b
    // %   x 
    // % 原问题
    // % Performs direct control allocation by solving the LP
    // %   max z=a   subj. to  Bu = av
    // %   a,u               umin <= u <= umax
    // % If a > 1, set u = u/a.
    // % Note: This function has not been optimized for speed.
    // %  Inputs:
    // %  -------
    // % B     control effectiveness matrix (k x m)
    // % v     commanded virtual control (k x 1)
    // % umin  lower position limits (m x 1)
    // % umax  upper position limits (m x 1)
    // %  Outputs:
    // %  -------
    // % u     optimal control (m x 1)
    // % a     scaling factor  
    // %% 整理成
    // %   min z=[0 -1]x   subj. to  [B -v]x = 0
    // %   x                       [I; -I]x <= [umax; Inf; -umin; 0] 
    // %   其中 x=[u; a]
    // % 对应《凸优化》p139,记为
    // %   min z=c*x   subj. to  Aeq*x = beq
    // %   x                     G*x <= h
    // But for using ALGLIB, the statement is:
    // %   min z=[0 -1]x   
    // %   subj. to  [B -v]x = 0  
    //             [umin; 0] <=x <= [umax; Inf] 
    // Variable is x=[u; a].
    // If a > 1, set u = u/a, 
    // B=[-0.5     0       0.5     0;
    //     0      -0.5     0       0.5;
    //     0.25    0.25    0.25    0.25];
        real_2d_array A = "[[-0.5, 0, 0.5, 0, 0.2],[0, -0.5, 0, 0.5, 0.1], [0.25, 0.25, 0.25, 0.25, 0]]"; 
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


        // Solve again
        real_1d_array minus_v = "[0.2, 0.1, 0.1]"; //v= -[0.2 0.1 0.1];
        A[0][4]=minus_v[0];
        A[1][4]=minus_v[1];
        A[2][4]=minus_v[2];
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
        double y_all1[3]={(double) -0.2, (double) -0.1, (double)  -0.0};
        double y_all[3]={(double) -0.2, (double) -0.1, (double)  -0.1};
        double u_all[4]={(double) 0.0, (double) 0.0, (double)  0.0, (double)  0.0};;
        double z_all=(double) 0.0;
        double iters_all=(double) 0.0;
        double _uMin[4] ={};
	    double _uMax[4] ={};
        for (int i = 0; i < 4; i++)
		{
			_uMin[i] = (double) -0.3491;
			_uMax[i] = (double) 0.3491;
		}
        dir_alloc_sim(y_all1, _uMin, _uMax, u_all, &z_all, &iters_all);
        start = std::chrono::high_resolution_clock::now();
        dir_alloc_sim(y_all, _uMin, _uMax, u_all, &z_all, &iters_all);
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        std::cout << "dir_alloc_sim execution time: " << elapsed.count() << "s\n";
        // 使用循环打印数组元素
        std::cout << "u_all: " << " ";
        for (int i = 0; i < 4; ++i) {
            std::cout << u_all[i] << " ";
        }
        std::cout << std::endl; 


    // }
    // catch(alglib::ap_error alglib_exception)
    // {
    //     printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
    //     return 1;
    // }
    return 0;
}
