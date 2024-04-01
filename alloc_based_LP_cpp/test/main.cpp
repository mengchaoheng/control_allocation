#include <iostream>
#include <matrix/math.hpp>
#include"ControlAllocation.h"

using namespace matrix;


int main() {
    
    
    float _B[3][4] = { {-0.4440,0.0,0.4440,0.0}, {0.0,-0.4440,0.0,0.4440},{0.2070,0.2070,0.2070,0.2070}};
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
    float yd[3]={0.2, -0.1, 0.1};
    

    // 线性规划数据
    LinearProgrammingProblem<3, 5> problem;
    problem.itlim = 100;
    //填数据
    problem.inB[0]=0;
    problem.inB[1]=1;
    problem.inB[2]=3;

    problem.e[0] = true;
    problem.e[1] = true;
    problem.e[2] = false;
    problem.e[3] = true;
    problem.e[4] = true;
    
    float upper_lam=1e4;
    for(int i=0; i<problem.m; ++i)
    {
        float temp=0;
        for(int j=0; j<problem.n-1; ++j)
        {
            problem.A[i][j] =_B[i][j];
            temp +=-_B[i][j]*_uMin[j];
        }
        problem.A[i][problem.n-1] =-yd[i];
        problem.b[i] = temp;
    }
    for(int i=0; i<problem.n-1; ++i)
    {
        problem.c[i] =0;
    }
    problem.c[problem.n-1] =-1;
    for(int i=0; i<problem.n; ++i)
    {
        problem.h[i] =_uMax[i]-_uMin[i];
    }
    problem.h[problem.n-1]=upper_lam;

    // 调用函数模板
    auto start = std::chrono::high_resolution_clock::now();
    LinearProgrammingResult<3, 5> result = BoundedRevisedSimplex(problem);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "execution time: " << elapsed.count() << "s\n";
    
    // 使用结果
    // result.y0, result.inB, result.e, result.errout
    int err = 0;
    float xout[problem.n];
    for(int i=0;i<problem.n;++i){
        xout[i]=0;
    }
    for(int i=0;i<problem.m;++i){
        xout[result.inB[i]]=result.y0[i];
    }
    for(int i=0;i<problem.n;++i){
        if(!result.e[i]){
            xout[i]=-xout[i]+problem.h[i];
        }
    }
    if(result.itlim<=0){
        err = 3;
        std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
    }
	if(result.errout)
    {
        err = 1;
        std::cout << "Solver error"<< std::endl;
    }
    float u_px4_matrix[problem.n-1];
    for(int i=0;i<problem.n-1;++i){
        u_px4_matrix[i]=xout[i]+_uMin[i];
    }
    if(xout[problem.n-1]>1){
        for(int i=0;i<problem.n-1;++i){
            u_px4_matrix[i]/=xout[problem.n-1];
        }
    }
    
    std::cout << "u_px4_matrix: [";
    for (size_t i = 0; i < problem.n-1; ++i) {
        std::cout << u_px4_matrix[i];
        if (i < problem.n-2) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    //飞机数据
    // 示例代码
    // 使用模板类时，可以指定 ControlSize 和 EffectSize 的具体值
    // 例如：
    float l1=0.148;float l2=0.069;float k_v=3;
    Aircraft<3, 4> df_4(_B, -0.3491, 0.3491); // 创建一个具有 4 个操纵向量和 3 个广义力矩的飞行器对象
    // 分配器数据：
    DP_LP_ControlAllocator<3, 4> DP_LPCA(df_4, yd); // 创建一个控制分配器对象，用于具有 4 个操纵向量和 3 个广义力矩的飞行器(转化为线性规划问题，其维数和参数 <3, 4> 有关。)
    // 然后可以使用飞行器对象和控制分配器对象进行操作
    // 调用 u = DP_LP_ControlAllocator.allocateControl(yd)
    float* u = new float[4];
    start = std::chrono::high_resolution_clock::now();
    u = DP_LPCA.allocateControl(yd);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    std::cout << "DP_LPCA.allocateControl execution time: " << elapsed.count() << "s\n";
    std::cout << "u: [";
    for (size_t i = 0; i < 4; ++i) {
        std::cout << u[i];
        if (i < 3) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    return 0;
}
