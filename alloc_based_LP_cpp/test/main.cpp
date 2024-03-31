#include <iostream>
#include <matrix/math.hpp>
using namespace matrix;
// 定义线性规划问题结构体
template<int M, int N>
struct LinearProgrammingProblem {
    
    int m=M;
    int n=N;
    int inB[M];
    int inD[N-M];
    int itlim;
    float A[M][N];
    float b[N];
    float c[N];
    float h[N];
    bool e[N];
    // 函数用于计算两个正整数集合的差
    void setdiff(int setA[], int sizeA, int setB[], int sizeB, int result[]) {
        int sizeResult = 0;
        for (int i = 0; i < sizeA; ++i) {
            bool foundInB = false;
            // 检查当前 setA 中的元素是否在 setB 中
            for (int j = 0; j < sizeB; ++j) {
                if (setA[i] == setB[j]) {
                    foundInB = true;
                    break;
                }
            }
            // 如果当前元素不在 setB 中，则将其添加到结果中
            if (!foundInB) {
                result[sizeResult++] = setA[i];
            }
        }
    }
    int* generateSequence(int i, int n) {
        int* result = new int[n - i + 1]; // 动态分配数组内存

        for (int num = i, index = 0; num <= n; ++num, ++index) {
            result[index] = num;
        }

        return result;
    }
    
};

// 定义结果结构体
template<int M, int N>
struct LinearProgrammingResult {
    float y0[M];
    int inB[M];
    bool e[N];
    int itlim;
    bool errout;
    // 其他结果成员
};



// 定义函数模板
template<int M, int N>
LinearProgrammingResult<M, N> BoundedRevisedSimplex(LinearProgrammingProblem<M, N>& problem) {
    LinearProgrammingResult<M, N> result;
    // 实现线性规划算法
    // 使用 problem.inB, problem.inD, problem.itlim, problem.A, problem.b, problem.c, problem.h, problem.e
    float tol=1e-7;
    const int n_m=N-M;
    int* nind = problem.generateSequence(0, n_m-1);
    
    int* ind_all = problem.generateSequence(0, N-1);
    problem.setdiff(ind_all, N, problem.inB, M, problem.inD);

    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            if(!problem.e[j])
            {
                problem.A[i][j] *=-1;
                problem.b[i]+=problem.A[i][j]*problem.h[j];
            }
        }
    }
    for(int j=0; j<N; ++j)
    {
        if(!problem.e[j])
        {
            problem.c[j] *=-1;
        }
    }

    //==============================
    matrix::SquareMatrix<float, M> A_inB;
    matrix::Matrix<float, M, n_m> A_inD;
    matrix::Vector<float, M> c_inB;
    matrix::Vector<float, n_m> c_inD;
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<M; ++j)
        {
            A_inB(i,j)=problem.A[i][problem.inB[j]];
            if(j<n_m)
            {
                A_inD(i,j)=problem.A[i][problem.inD[j]];
            }
        }
        c_inB(i)=problem.c[problem.inB[i]];
    }
    for(int i=0; i<n_m; ++i)
    {
        c_inD(i)=problem.c[problem.inD[i]];
    }
    matrix::Vector<float, M> b_vec(problem.b);

    // inital some value
    Matrix<float, 1UL, M> lamt;
    lamt.setZero();
    Matrix<float, 1UL, n_m> rdt;
    rdt.setZero();
    matrix::Vector<float, M> A_qel;
    A_qel.setZero();
    matrix::Vector<float, M> yq;
    yq.setZero();
    matrix::Vector<float, M> rat;
    rat.setZero();
    
    //  %Initial Solution
    matrix::Vector<float, M> y0 = inv(A_inB)*b_vec;
    bool done = false;
    bool unbounded = false;
     while ((!done  || !unbounded ) && (problem.itlim > 0))
    {
        problem.itlim = problem.itlim-1;
        lamt= (inv(A_inB).transpose()*c_inB).transpose();
        rdt = c_inD.transpose()-lamt*A_inD;
        float minr;
        size_t qind;
        min(rdt.transpose(), minr, qind);
        if(minr >=0)  // If all relative costs are positive then the solution is optimal
        { 
            done = true;
            break;
        }
        int qel = problem.inD[qind];  // Unknown to Enter the basis minimizes relative cost
        A_qel(0)=problem.A[0][qel];
        A_qel(1)=problem.A[1][qel];
        A_qel(2)=problem.A[2][qel];
        yq=inv(A_inB)* A_qel; // Vector to enter in terms of the current Basis vector
        bool flag=false;
        
        for(int i=0;i<M;++i){
            if(std::abs(yq(i)) > tol)
            {
                flag = true; // Check this condition
                break;
            }
        }
        if(!flag)
        {
            unbounded = true; // Check this condition
            break;
        }
        // Recompute rations and determine variable to leave
        
        float hinB[M];
        for(int i=0;i<M;++i)
        {
            if(std::abs(yq(i))>tol)
            {
                rat(i)=y0(i)/yq(i);
                
            }
            else
            {
                rat(i)=INFINITY;
                /* code */
            }
        }
        for(int i=0;i<M;++i)
        {
            hinB[i]=problem.h[problem.inB[i]];
            if(yq(i)<0)
            {                    
                rat(i)-=hinB[i]/yq(i);
            }
        }
         // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        float minrat=rat(0);
        size_t p=0;
        min(rat, minrat, p);
        // If the minimum ratio is zero, then the solution is degenerate and the entering
        // variable will not change the basis---invoke Bland's selection rule to avoid
        // cycling.
        if (std::abs(minrat) <= tol)
        {
            //Find negative relative cost
            for(int i=0;i<n_m;++i)
            {
                if(rdt(1,i)<0){ //Note that since minr <0 indm is not empty   
                    qind=nind[i];
                    qel = problem.inD[qind];//Unknown to Enter the basis is first indexed to avoid cycling
                    break;
                }
            }
            A_qel(0)=problem.A[0][qel];
            A_qel(1)=problem.A[1][qel];
            A_qel(2)=problem.A[2][qel];
            yq=inv(A_inB)* A_qel;
            bool flag=false;
            for(int i=0;i<M;++i){
                if(std::abs(yq(i)) > tol)
                {
                    flag = true; // Check this condition
                    break;
                }
            }
            if(!flag)
            {
                unbounded = true; // Check this condition
                // break;
            }
            // Recompute rations and determine variable to leave
            // Recompute rations and determine variable to leave
            float hinB[M];
            for(int i=0;i<M;++i)
            {
                hinB[i]=problem.h[problem.inB[i]];
                if(std::abs(yq(i))>tol)
                {
                    rat(i)=y0(i)/yq(i);
                    if(yq(i)<0)
                    {                    
                        rat(i)-=hinB[i]/yq(i);
                    }
                }
                else
                {
                    rat(i)=INFINITY;
                    /* code */
                }
            }
            // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
            minrat=rat(0);
            p=0;
            min(rat, minrat, p);
        }
        if (minrat >= problem.h[qel])
        {
            // std::cout << " Case 1 "<< std::endl; 
            problem.e[qel] =!problem.e[qel];
            for(int i=0; i<M; ++i)
            {
                problem.A[i][qel] *= -1;
                b_vec(i)+=problem.A[i][qel]*problem.h[qel];
            }
            problem.c[qel] *= -1;

        }
        else if(yq(p) > 0)
        {
            // std::cout << " Case 21 "<< std::endl; 
            int pel = problem.inB[p];
            problem.inB[p]= qel;
            problem.inD[qind]= pel;
            // update x_inX
            for(int i=0; i<M; ++i)
            {
                A_inB(i,p)=problem.A[i][qel];
            }
            for(int i=0; i<n_m; ++i)
            {
                c_inB(p)=problem.c[qel];
            }
            for(int i=0; i<M; ++i)
            {
                A_inD(i,qind)=problem.A[i][pel];
            }
            for(int i=0; i<n_m; ++i)
            {
                c_inD(qind)=problem.c[pel];
            }
        }
        else
        {
            // std::cout << " Case 22 "<< std::endl; 
            int pel = problem.inB[p];
            problem.e[pel]=!problem.e[pel];
            for(int i=0; i<M; ++i)
            {
                problem.A[i][pel] *= -1;
                b_vec(i)+=problem.A[i][pel]*problem.h[pel];
            }
            problem.inB[p]= qel;
            problem.inD[qind]= pel;
            problem.c[pel] *= -1;
            // update x_inX
            for(int i=0; i<M; ++i)
            {
                A_inB(i,p)=problem.A[i][qel];
            }
            
            for(int i=0; i<n_m; ++i)
            {
                c_inB(p)=problem.c[qel];
            }
            for(int i=0; i<M; ++i)
            {
                A_inD(i,qind)=problem.A[i][pel];
            }
            for(int i=0; i<n_m; ++i)
            {
                c_inD(qind)=problem.c[pel];
            }
        }
        y0=inv(A_inB)* b_vec;
    }
    result.errout = unbounded; 
    // 设置 result.y0, result.inB, result.e 等结果
    for(int i=0; i<M; ++i)
    {
        result.y0[i]=y0(i);
        result.inB[i]=problem.inB[i];
    }
    for(int i=0; i<N; ++i)
    {
        result.e[i]=problem.e[i];
    }
    result.itlim=problem.itlim;
    return result;
}


int main() {
    //飞机数据
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
    // 分配器数据：

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
    float yd[3]={0.2, -0.1, 0.1};
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
    LinearProgrammingResult<3, 5>result = BoundedRevisedSimplex(problem);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << " execution time: " << elapsed.count() << "s\n";
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
    return 0;
}
