
#include <matrix/math.hpp>
#include <iostream>

using namespace matrix;

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
// 定义线性规划问题结构体
template<int M, int N>
struct LinearProgrammingProblem {
    int m=M;
    int n=N;
    int inB[M];
    int inD[N-M]; // 
    int itlim;
    float A[M][N];
    float b[N];
    float c[N];
    float h[N];
    bool e[N];
    float tol=1e-8;
    // 默认构造函数，将所有成员变量初始化为0
    LinearProgrammingProblem() : m(M), n(N), itlim(0) {
        // 将数组成员变量初始化为0
        for (int i = 0; i < M; ++i) {
            inB[i] = 0;
        }
        for (int i = 0; i < N-M; ++i) {
            inD[i] = 0;
        }
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                A[i][j] = 0.0f;
            }
        }
        for (int i = 0; i < N; ++i) {
            b[i] = 0.0f;
            c[i] = 0.0f;
            h[i] = 0.0f;
            e[i] = false;
        }
    }
    // 拷贝构造函数
    LinearProgrammingProblem(const LinearProgrammingProblem<M, N>& other) {
        m = other.m;
        n = other.n;
        itlim = other.itlim;
        tol = other.tol;
        // 复制数组成员变量
        for (int i = 0; i < M; ++i) {
            inB[i] = other.inB[i];
        }
        for (int i = 0; i < N - M; ++i) {
            inD[i] = other.inD[i];
        }
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                A[i][j] = other.A[i][j];
            }
        }
        for (int i = 0; i < N; ++i) {
            b[i] = other.b[i];
            c[i] = other.c[i];
            h[i] = other.h[i];
            e[i] = other.e[i];
        }
    }
};
// 定义结果结构体
template<int M, int N>
struct LinearProgrammingResult {
    float y0[M];
    int inB[M];
    bool e[N];
    int iters;
    bool errout;
    // 其他结果成员
    // 默认构造函数，将所有成员变量初始化为0
    LinearProgrammingResult() : iters(0), errout(false) {
        // 将数组成员变量初始化为0
        for (int i = 0; i < M; ++i) {
            y0[i] = 0.0f;
            inB[i] = 0;
        }
        for (int i = 0; i < N; ++i) {
            e[i] = false;
        }
    }
    // 拷贝构造函数
    LinearProgrammingResult(const LinearProgrammingResult<M, N>& other) {
        // 将成员变量从另一个对象复制到当前对象
        for (int i = 0; i < M; ++i) {
            y0[i] = other.y0[i];
            inB[i] = other.inB[i];
        }
        for (int i = 0; i < N; ++i) {
            e[i] = other.e[i];
        }
        iters = other.iters;
        errout = other.errout;
    }
};

template<int M, int N>
class LinearProgrammingSolverForAC {
public:
    // 默认构造函数，初始化所有成员变量
    LinearProgrammingSolverForAC() 
        : problem(),                // 使用默认构造函数初始化 problem
          result(),                 // 使用默认构造函数初始化 result
          tol(1e-7),                // 将 tol 初始化为 1e-7
          n_m(N - M)               // 计算并初始化 n_m
    {
        // nind = generateSequence(0, N-M-1);
        // ind_all = generateSequence(0, N-1);
        for (int i = 0; i < N - M; ++i) {
            nind[i] = i; // 使用 i 初始化 nind 数组的元素
        }
        // ind_all = generateSequence(0, N-1);
        // 初始化 ind_all 数组
        for (int i = 0; i < N ; ++i) {
            ind_all[i] = i; // 使用 i 初始化 nind 数组的元素
        }
        A_inB.setZero();
        A_inD.setZero();
        c_inB.setZero();
        c_inD.setZero();
        lamt.setZero();
        rdt.setZero();
        A_qel.setZero();
        yq.setZero();
        rat.setZero();
    }
    // 构造函数，接受一个 LinearProgrammingProblem 对象作为参数
    LinearProgrammingSolverForAC(const LinearProgrammingProblem<M, N> inputProblem) : problem(inputProblem) {
        // nind = generateSequence(0, N-M-1);
        // 初始化 nind 数组
        for (int i = 0; i < N - M; ++i) {
            nind[i] = i; // 使用 i 初始化 nind 数组的元素
        }
        // ind_all = generateSequence(0, N-1);
        // 初始化 ind_all 数组
        for (int i = 0; i < N ; ++i) {
            ind_all[i] = i; // 使用 i 初始化 nind 数组的元素
        }
        A_inB.setZero();
        A_inD.setZero();
        c_inB.setZero();
        c_inD.setZero();
        lamt.setZero();
        rdt.setZero();
        A_qel.setZero();
        yq.setZero();
        rat.setZero();
    }
    // 成员函数调用 Simplex （改自 BoundedRevisedSimplex ），并返回结果
    // void solve() {
    //     result = Simplex(problem);
    // }
    LinearProgrammingResult<M, N>  solve() {
        result = Simplex(problem);
        return result;
    }
    LinearProgrammingProblem<M, N> problem;
    LinearProgrammingResult<M, N> result;
    // for simplex
    //==============================
    matrix::SquareMatrix<float, M> A_inB;
    matrix::Matrix<float, M, N-M> A_inD;
    matrix::Vector<float, M> c_inB;
    matrix::Vector<float, N-M> c_inD;
    // inital some value
    Matrix<float, 1UL, M> lamt;
    Matrix<float, 1UL, N-M> rdt;
    matrix::Vector<float, M> A_qel;
    matrix::Vector<float, M> yq;
    matrix::Vector<float, M> rat;
    float tol=1e-7;
    const int n_m=N-M;
    int nind[N-M]; 
    int ind_all[N]; 
    LinearProgrammingResult<M, N> Simplex(LinearProgrammingProblem<M, N> problem){
        LinearProgrammingResult<M, N> result;
        // 实现线性规划算法
        // 使用 problem.inB, problem.inD, problem.itlim, problem.A, problem.b, problem.c, problem.h, problem.e
        // Partition A
        std::cout << "ind_all: [";
        for (size_t i = 0; i <N; ++i) {
            std::cout << ind_all[i];
            if (i <N- 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
        setdiff(ind_all, N, problem.inB, M, problem.inD);
        std::cout << "problem.inB: [";
        for (size_t i = 0; i <M; ++i) {
            std::cout << problem.inB[i];
            if (i <M- 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
        std::cout << "problem.inD: [";
        for (size_t i = 0; i <N-M; ++i) {
            std::cout << problem.inD[i];
            if (i <N-M- 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;


        std::cout << "A:" << std::endl;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                std::cout << problem.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "b: [";
        for (size_t i = 0; i < M; ++i) {
            std::cout << problem.b[i];
            if (i < M - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;

        std::cout << "c: [";
        for (size_t i = 0; i < N; ++i) {
            std::cout << problem.c[i];
            if (i < N - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;

        std::cout << "h: [";
        for (size_t i = 0; i < N; ++i) {
            std::cout << problem.h[i];
            if (i < N- 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
        // Adjust signs problem if variables are initialized at upper bounds.
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
        std::cout << "A:" << std::endl;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                std::cout << problem.A[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "b: [";
        for (size_t i = 0; i < M; ++i) {
            std::cout << problem.b[i];
            if (i < M - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;

        std::cout << "c: [";
        for (size_t i = 0; i < N; ++i) {
            std::cout << problem.c[i];
            if (i < N - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;

        std::cout << "h: [";
        for (size_t i = 0; i < N; ++i) {
            std::cout << problem.h[i];
            if (i < N- 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
        //==============================

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
        // Initial Solution
        // matrix::Vector<float, M> y0 = inv(A_inB)*b_vec;
        matrix::LeastSquaresSolver<float, M,M> LSsolver0(A_inB);
        matrix::Vector<float, M> y0 = LSsolver0.solve(b_vec);
        bool done = false;
        bool unbounded = false;
        int iters =0;
        while ((!done  || !unbounded ) && (iters <= problem.itlim))
        {
            iters = iters+1;
            // lamt= (inv(A_inB).transpose()*c_inB).transpose();
            matrix::LeastSquaresSolver<float, M,M> LSsolver_lamt(A_inB.transpose());
            lamt = LSsolver_lamt.solve(c_inB).transpose();
            rdt = c_inD.transpose()-lamt*A_inD;
            std::cout << "lamt:";
            lamt.print();
            std::cout << "rdt:";
            rdt.print();
            float minr;
            size_t qind;
            min(rdt.transpose(), minr, qind);
            std::cout << "minr:"<<minr<<std::endl;
            std::cout << "qind:"<<qind<<std::endl;
            if(minr >=0)  // If all relative costs are positive then the solution is optimal
            { 
                done = true;
                break;
            }
            int qel = problem.inD[qind];  // Unknown to Enter the basis minimizes relative cost
            A_qel(0)=problem.A[0][qel];
            A_qel(1)=problem.A[1][qel];
            A_qel(2)=problem.A[2][qel];
            std::cout << "qel:"<<qel<<std::endl;

            // yq=inv(A_inB)* A_qel; // Vector to enter in terms of the current Basis vector
            matrix::LeastSquaresSolver<float, M,M> LSsolver1(A_inB);
            yq = LSsolver1.solve(A_qel);
            std::cout << "yq:";
            yq.print();
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
                if(yq(i)<0 && std::abs(yq(i))>tol )
                {                    
                    rat(i)-=hinB[i]/yq(i);
                }
            }
            std::cout << "rat:";
            rat.print();
            // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
            float minrat=rat(0);
            size_t p=0;
            min(rat, minrat, p);
            std::cout << "minrat:"<<minrat<<std::endl;
            std::cout << "p:"<<p<<std::endl;
            // If the minimum ratio is zero, then the solution is degenerate and the entering
            // variable will not change the basis---invoke Bland's selection rule to avoid
            // cycling.
            if (std::abs(minrat) <= tol)
            {
                //Find negative relative cost
                std::cout << " Find negative relative cost "<< std::endl; 
                for(int i=0;i<N-M;++i)
                {
                    // std::cout << "rdt(0,i):"<<rdt(0,i)<<std::endl; 
                    if(rdt(0,i)<0){ //Note that since minr <0 indm is not empty 
                        qind=nind[i];
                        qel = problem.inD[qind];//Unknown to Enter the basis is first indexed to avoid cycling
                        std::cout << "qind:"<<qind<<std::endl;
                        std::cout << "qel:"<<qel<<std::endl;
                        break;
                    }
                }
                A_qel(0)=problem.A[0][qel];
                A_qel(1)=problem.A[1][qel];
                A_qel(2)=problem.A[2][qel];
                // yq=inv(A_inB)* A_qel;
                matrix::LeastSquaresSolver<float, M,M> LSsolver2(A_inB);
                yq = LSsolver2.solve(A_qel);
                std::cout << "yq:";
                yq.print();
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
                std::cout << "rat:";
                rat.print();
                // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
                minrat=rat(0);
                p=0;
                min(rat, minrat, p);
                std::cout << "minrat:"<<minrat<<std::endl;
                std::cout << "p:"<<p<<std::endl;
            }
            if (minrat >= problem.h[qel])
            {
                std::cout << " Case 1 "<< std::endl; 
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
                std::cout << " Case 21 "<< std::endl; 
                int pel = problem.inB[p];
                std::cout << "pel:"<<pel<<std::endl;
                std::cout << "qel:"<<qel<<std::endl;
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
                std::cout << "A:" << std::endl;
                for (size_t i = 0; i < M; ++i) {
                    for (size_t j = 0; j < N; ++j) {
                        std::cout << problem.A[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "b: [";
                for (size_t i = 0; i < M; ++i) {
                    std::cout << problem.b[i];
                    if (i < M - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
                std::cout << "b_vec:";
                b_vec.print();
                std::cout << "problem.inB: [";
                for (size_t i = 0; i <M; ++i) {
                    std::cout << problem.inB[i];
                    if (i <M- 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
                std::cout << "problem.inD: [";
                for (size_t i = 0; i <N-M; ++i) {
                    std::cout << problem.inD[i];
                    if (i <N-M- 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
            }
            else
            {
                std::cout << " Case 22 "<< std::endl; 
                int pel = problem.inB[p];
                std::cout << "pel:"<<pel<<std::endl;
                problem.e[pel]=!problem.e[pel];
                std::cout << "problem.e[pel]:"<<problem.e[pel]<<std::endl;
                for(int i=0; i<M; ++i)
                {
                    problem.A[i][pel] *= -1;
                    b_vec(i)+=problem.A[i][pel]*problem.h[pel];
                    problem.b[i]=b_vec(i);

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
                std::cout << "A:" << std::endl;
                for (size_t i = 0; i < M; ++i) {
                    for (size_t j = 0; j < N; ++j) {
                        std::cout << problem.A[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "b: [";
                for (size_t i = 0; i < M; ++i) {
                    std::cout << problem.b[i];
                    if (i < M - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
                std::cout << "b_vec:";
                b_vec.print();
                std::cout << "c: [";
                for (size_t i = 0; i < N; ++i) {
                    std::cout << problem.c[i];
                    if (i < N - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;

                std::cout << "problem.inB: [";
                for (size_t i = 0; i <M; ++i) {
                    std::cout << problem.inB[i];
                    if (i <M- 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
                std::cout << "problem.inD: [";
                for (size_t i = 0; i <N-M; ++i) {
                    std::cout << problem.inD[i];
                    if (i <N-M- 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << "]" << std::endl;
            }
            // y0=inv(A_inB)* b_vec;

            matrix::LeastSquaresSolver<float, M,M> LSsolver(A_inB);
            y0 = LSsolver.solve(b_vec);
            std::cout << "y0:";
            y0.print();
            std::cout << "A_inB:";
            A_inB.print();
        }
        result.errout = unbounded; 
        // y0.print();
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
        result.iters=iters;
        std::cout << "iters:"<<iters<<std::endl;
        std::cout << "problem.itlim:"<<problem.itlim<<std::endl;
        return result;
    }
    // 获取结果
    const LinearProgrammingResult<M, N>& getResult() const {
        return result;
    }
};


// 定义函数模板
template<int M, int N>
LinearProgrammingResult<M, N> BoundedRevisedSimplex(LinearProgrammingProblem<M, N> problem) {
    LinearProgrammingResult<M, N> result;
    // 实现线性规划算法
    // 使用 problem.inB, problem.inD, problem.itlim, problem.A, problem.b, problem.c, problem.h, problem.e
    // float tol=1e-7;
    const int n_m=N-M;
    int* nind = generateSequence(0, n_m-1);
    // std::cout << "nind: [";
    // for (size_t i = 0; i <N- M; ++i) {
    //     std::cout << nind[i];
    //     if (i <N- M - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    //=============
    int* ind_all = generateSequence(0, N-1);
    //=============
    // std::cout << "ind_all: [";
    // for (size_t i = 0; i <N; ++i) {
    //     std::cout << ind_all[i];
    //     if (i <N- 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    //=============
    setdiff(ind_all, N, problem.inB, M, problem.inD);
    //=============
    // std::cout << "problem.inB: [";
    // for (size_t i = 0; i <M; ++i) {
    //     std::cout << problem.inB[i];
    //     if (i <M- 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;
    // std::cout << "problem.inD: [";
    // for (size_t i = 0; i <N-M; ++i) {
    //     std::cout << problem.inD[i];
    //     if (i <N-M- 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;


    // std::cout << "A:" << std::endl;
    // for (size_t i = 0; i < M; ++i) {
    //     for (size_t j = 0; j < N; ++j) {
    //         std::cout << problem.A[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "b: [";
    // for (size_t i = 0; i < M; ++i) {
    //     std::cout << problem.b[i];
    //     if (i < M - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    // std::cout << "c: [";
    // for (size_t i = 0; i < N; ++i) {
    //     std::cout << problem.c[i];
    //     if (i < N - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    // std::cout << "h: [";
    // for (size_t i = 0; i < N; ++i) {
    //     std::cout << problem.h[i];
    //     if (i < N- 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    
    //

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
    // std::cout << "A:" << std::endl;
    // for (size_t i = 0; i < M; ++i) {
    //     for (size_t j = 0; j < N; ++j) {
    //         std::cout << problem.A[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "b: [";
    // for (size_t i = 0; i < M; ++i) {
    //     std::cout << problem.b[i];
    //     if (i < M - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    // std::cout << "c: [";
    // for (size_t i = 0; i < N; ++i) {
    //     std::cout << problem.c[i];
    //     if (i < N - 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    // std::cout << "h: [";
    // for (size_t i = 0; i < N; ++i) {
    //     std::cout << problem.h[i];
    //     if (i < N- 1) {
    //         std::cout << ", ";
    //     }
    // }
    // std::cout << "]" << std::endl;

    
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
    Matrix<float, 1, M> lamt;
    lamt.setZero();
    Matrix<float, 1, n_m> rdt;
    rdt.setZero();
    matrix::Vector<float, M> A_qel;
    A_qel.setZero();
    matrix::Vector<float, M> yq;
    yq.setZero();
    matrix::Vector<float, M> rat;
    rat.setZero();
    
    //  %Initial Solution
    // matrix::Vector<float, M> y0 = inv(A_inB)*b_vec;
    matrix::LeastSquaresSolver<float, M,M> LSsolver0(A_inB);
    matrix::Vector<float, M> y0 = LSsolver0.solve(b_vec);
    // std::cout << "y0:";
    // y0.print();
    bool done = false;
    bool unbounded = false;
    int iters =0;
     while ((!done  || !unbounded ) && (iters <= problem.itlim))
    {
        iters = iters+1;
        // lamt= (inv(A_inB).transpose()*c_inB).transpose();
        matrix::LeastSquaresSolver<float, M,M> LSsolver_lamt(A_inB.transpose());
        lamt = LSsolver_lamt.solve(c_inB).transpose();
        rdt = c_inD.transpose()-lamt*A_inD;
        // std::cout << "lamt:";
        // lamt.print();
        // std::cout << "rdt:";
        // rdt.print();
        float minr;
        size_t qind;
        min(rdt.transpose(), minr, qind);
        // std::cout << "minr:"<<minr<<std::endl;
        // std::cout << "qind:"<<qind<<std::endl;
        if(minr >=0)  // If all relative costs are positive then the solution is optimal. have to compare with 0 !
        { 
            done = true;
            break;
        }
        int qel = problem.inD[qind];  // Unknown to Enter the basis minimizes relative cost
        A_qel(0)=problem.A[0][qel];
        A_qel(1)=problem.A[1][qel];
        A_qel(2)=problem.A[2][qel];
        // std::cout << "qel:"<<qel<<std::endl;

        // yq=inv(A_inB)* A_qel; // Vector to enter in terms of the current Basis vector
        matrix::LeastSquaresSolver<float, M,M> LSsolver1(A_inB);
        yq = LSsolver1.solve(A_qel);
        // std::cout << "yq:";
        // yq.print();
        bool flag=false;
        
        for(int i=0;i<M;++i){
            if(std::abs(yq(i)) > problem.tol)
            {
                flag = true; // Check this condition
                break;
            }
        }
        if(!flag)
        {
            unbounded = true; // Check this condition
            std::cout << "simplex loop Solution is unbounded"<< std::endl; 
            break;
        }
        // Recompute rations and determine variable to leave
        // careful here
        float hinB[M];
        for(int i=0;i<M;++i)
        {
            if(std::abs(yq(i))>problem.tol)
            {
                rat(i)=y0(i)/yq(i);
                hinB[i]=problem.h[problem.inB[i]];
                if(yq(i)<0 ) // have to be compare with 0!!!
                {
                    // std::cout << "yq(i)<0 i:"<< i<< std::endl;                  
                    rat(i)-=hinB[i]/yq(i);
                }
            }
            else
            {
                rat(i)=INFINITY;
                /* code */
            }
        }
        // std::cout << "rat:";
        // rat.print();
         // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        float minrat=rat(0);
        size_t p=0;
        min(rat, minrat, p);
        // std::cout << "minrat:"<<minrat<<std::endl;
        // std::cout << "p:"<<p<<std::endl;
        // If the minimum ratio is zero, then the solution is degenerate and the entering
        // variable will not change the basis---invoke Bland's selection rule to avoid
        // cycling.
        if (std::abs(minrat) <= problem.tol)
        {
            //Find negative relative cost
            // std::cout << " Find negative relative cost "<< std::endl; 
            for(int i=0;i<N-M;++i)
            {
                // std::cout << "rdt(0,i):"<<rdt(0,i)<<std::endl; 
                if(rdt(0,i)<0){ //Note that since minr <0 indm is not empty 
                    qind=nind[i];
                    qel = problem.inD[qind];//Unknown to Enter the basis is first indexed to avoid cycling
                    // std::cout << "qind:"<<qind<<std::endl;
                    // std::cout << "qel:"<<qel<<std::endl;
                    break;
                }
            }
            A_qel(0)=problem.A[0][qel];
            A_qel(1)=problem.A[1][qel];
            A_qel(2)=problem.A[2][qel];
            // yq=inv(A_inB)* A_qel;
            matrix::LeastSquaresSolver<float, M,M> LSsolver2(A_inB);
            yq = LSsolver2.solve(A_qel);
            // std::cout << "yq:";
            // yq.print();
            bool flag=false;
            for(int i=0;i<M;++i){
                if(std::abs(yq(i)) > problem.tol)
                {
                    flag = true; // Check this condition
                    break;
                }
            }
            if(!flag)
            {
                unbounded = true; // Check this condition
                std::cout << "simplex loop Solution is unbounded"<< std::endl; 
                break;
            }
            // Recompute rations and determine variable to leave
            // Recompute rations and determine variable to leave
            float hinB[M];
            for(int i=0;i<M;++i)
            {
                hinB[i]=problem.h[problem.inB[i]];
                if(std::abs(yq(i))>problem.tol)
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
            // std::cout << "rat:";
            // rat.print();
            // Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
            minrat=rat(0);
            p=0;
            min(rat, minrat, p);
            // std::cout << "minrat:"<<minrat<<std::endl;
            // std::cout << "p:"<<p<<std::endl;
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

            for(int i=0; i<M; ++i)
            {
                A_inD(i,qind)=problem.A[i][qel];
            }
            c_inD(qind)=problem.c[qel];
       

        }
        else if(yq(p) > 0)
        {
            // std::cout << " Case 21 "<< std::endl; 
            int pel = problem.inB[p];
            // std::cout << "pel:"<<pel<<std::endl;
            // std::cout << "qel:"<<qel<<std::endl;
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
            // std::cout << "pel:"<<pel<<std::endl;
            problem.e[pel]=!problem.e[pel];
            // std::cout << "problem.e[pel]:"<<problem.e[pel]<<std::endl;
            for(int i=0; i<M; ++i)
            {
                problem.A[i][pel] *= -1;
                b_vec(i)+=problem.A[i][pel]*problem.h[pel];
                problem.b[i]=b_vec(i);

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
        // update x_inX move to here
        //     for(int i=0; i<M; ++i)
        // {
        //     for(int j=0; j<M; ++j)
        //     {
        //         A_inB(i,j)=problem.A[i][problem.inB[j]];
        //         if(j<n_m)
        //         {
        //             A_inD(i,j)=problem.A[i][problem.inD[j]];
        //         }
        //     }
        //     c_inB(i)=problem.c[problem.inB[i]];
        // }
        // for(int i=0; i<n_m; ++i)
        // {
        //     c_inD(i)=problem.c[problem.inD[i]];
        // }
        // y0=inv(A_inB)* b_vec;

        matrix::LeastSquaresSolver<float, M,M> LSsolver(A_inB);
        y0 = LSsolver.solve(b_vec);
        // std::cout << "y0:";
        // y0.print();
        // std::cout << "A_inB:";
        // A_inB.print();
        // std::cout << "A_inD:";
        // A_inD.print();
        // std::cout << "c_inB:";
        // c_inB.print();
        // std::cout << "c_inD:";
        // c_inD.print();
        
        // std::cout << "A:" << std::endl;
        // for (size_t i = 0; i < M; ++i) {
        //     for (size_t j = 0; j < N; ++j) {
        //         std::cout << problem.A[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "b: [";
        // for (size_t i = 0; i < M; ++i) {
        //     std::cout << problem.b[i];
        //     if (i < M - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "b_vec:";
        // b_vec.print();
        // std::cout << "c: [";
        // for (size_t i = 0; i < N; ++i) {
        //     std::cout << problem.c[i];
        //     if (i < N - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        // std::cout << "problem.inB: [";
        // for (size_t i = 0; i <M; ++i) {
        //     std::cout << problem.inB[i];
        //     if (i <M- 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "problem.inD: [";
        // for (size_t i = 0; i <N-M; ++i) {
        //     std::cout << problem.inD[i];
        //     if (i <N-M- 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;




        
    }
    result.errout = unbounded; 
    // y0.print();
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
    result.iters=iters;
    // std::cout << "iters:"<<iters<<std::endl;
    // std::cout << "problem.itlim:"<<problem.itlim<<std::endl;
    return result;
}




// 飞行器基类模板
template <int ControlSize, int EffectorSize>
class AircraftBase {
public:
    float controlVector[EffectorSize]; // 操纵向量
    float controlEffectMatrix[ControlSize][EffectorSize]; // 控制效应矩阵 (generalizedMomentSize X controlVectorSize)
    float upperLimits[EffectorSize]; // 操纵向量上限变量
    float lowerLimits[EffectorSize]; // 操纵向量下限变量
    // 构造函数
    // 拷贝构造函数
    AircraftBase(const AircraftBase& other) {
        // 复制 controlVector
        for (int i = 0; i < EffectorSize; ++i) {
            controlVector[i] = other.controlVector[i];
            upperLimits[i] = other.upperLimits[i];
            lowerLimits[i] = other.lowerLimits[i];
        }

        // 复制 controlEffectMatrix
        for (int i = 0; i < ControlSize; ++i) {
            for (int j = 0; j < EffectorSize; ++j) {
                controlEffectMatrix[i][j] = other.controlEffectMatrix[i][j];
            }
        }
    }
    AircraftBase() {
        // 初始化 controlVector, controlEffectMatrix, generalizedMoment, upperLimits, lowerLimits 等数组
        // 可以使用默认初始化或者自定义初始化方式
        // 例如：
        for (int i = 0; i < EffectorSize; ++i) {
            controlVector[i] = 0.0f;
            upperLimits[i] = 0.0f;
            lowerLimits[i] = 0.0f;
            for (int j = 0; j < ControlSize; ++j) {
                controlEffectMatrix[j][i] = 0.0f;
            }
        }
    }
    // 析构函数
    ~AircraftBase() {
        // 不需要手动释放内存，因为数组是在栈上分配的，会在对象生命周期结束时自动释放
    }
};

// 飞行器类模板, 不同飞机定义新的类继承基类
template <int ControlSize, int EffectorSize>
class Aircraft : public AircraftBase<ControlSize, EffectorSize> {
private:
    // 添加特定飞行器类型的模型参数
public:
    // 构造函数
    // 拷贝构造函数
    // 拷贝构造函数
    Aircraft(const Aircraft& other) : AircraftBase<ControlSize, EffectorSize>(other) {
        // 将其他对象的成员变量值复制到新对象中
        l1 = other.l1;
        l2 = other.l2;
        k_v = other.k_v;
        upper = other.upper;
        lower = other.lower;
    }
    Aircraft() : AircraftBase<ControlSize, EffectorSize>() {
        // 可选的初始化代码
        l1=0;
        l2=0;
        k_v=0;
    }
    // 构造函数，接受对应于飞行器类模板参数的初始化参数
    Aircraft(float controlEffectMatrixInit[ControlSize][EffectorSize], 
             float upperLimitsInit[EffectorSize], 
             float lowerLimitsInit[EffectorSize]) {
        // 使用传入的初始化参数对飞行器的数组成员进行初始化
        for (int i = 0; i < EffectorSize; ++i) {
            this->controlVector[i] = 0;
            this->upperLimits[i] = upperLimitsInit[i];
            this->lowerLimits[i] = lowerLimitsInit[i];
            for (int j = 0; j < ControlSize; ++j) {
                this->controlEffectMatrix[j][i] = controlEffectMatrixInit[j][i];
            }
        }
    }
    // 构造函数，接受对应于飞行器类模板参数的初始化参数
    Aircraft(float upperLimitsInit[EffectorSize], 
             float lowerLimitsInit[EffectorSize]) {
        // 使用传入的初始化参数和飞行器
        for (int i = 0; i < EffectorSize; ++i) {
            this->controlVector[i] = 0;
            this->upperLimits[i] = upperLimitsInit[i];
            this->lowerLimits[i] = lowerLimitsInit[i];
            for (int j = 0; j < ControlSize; ++j) {
                this->controlEffectMatrix[j][i] = 0;
            }
        }
        // and define other value manual to set B (controlEffectMatrix).
    }
    Aircraft(float controlEffectMatrixInit[ControlSize][EffectorSize], float lower, float upper) : lower(lower), upper(upper) {
        // 使用传入的初始化参数和飞行器
        std::cout << "使用传入的controlEffectMatrixInit, lower, upper初始化参数和飞行器"<< std::endl; 
        for (int i = 0; i < EffectorSize; ++i) {
            this->controlVector[i] = 0;
            this->upperLimits[i] = upper;
            this->lowerLimits[i] = lower;
            for (int j = 0; j < ControlSize; ++j) {
                this->controlEffectMatrix[j][i] = controlEffectMatrixInit[j][i];
            }
        }
    }

    // 设置模型参数函数
    float l1;
    float l2;
    float k_v;
    float upper;
    float lower;

    // 析构函数
    ~Aircraft() {
        // 如果有需要释放的资源，可以在这里添加代码
    }

    // 其他成员函数和成员变量定义
};
// 控制分配基类模板
template <int ControlSize, int EffectorSize>
class ControlAllocatorBase {
public:
    // 构造函数
    ControlAllocatorBase() : aircraft() {
        // 在此初始化成员变量，或者留空
    }
    // 参数列表构造函数
    ControlAllocatorBase(const Aircraft<ControlSize, EffectorSize>& aircraft) 
        : aircraft(aircraft) { 
        // 使用传入的aircraft对象初始化aircraft成员
        std::cout << "ControlAllocatorBase使用传入的aircraft对象初始化aircraft成员"<< std::endl;
    }

    virtual float*  allocateControl(float input[ControlSize]) = 0;

    // 其他数学函数和成员变量定义
    Aircraft<ControlSize, EffectorSize> aircraft; // 构造函数设置
    
};
// 控制分配类模板
template <int ControlSize, int EffectorSize>
class DP_LP_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectorSize> {
private:
    // 添加算法设置参数
public:
    // 构造函数, 利用aircraft 预设置LinearProgrammingProblem, 再初始化 LinearProgrammingSolverForAC
    // 构造函数
    DP_LP_ControlAllocator(const Aircraft<ControlSize, EffectorSize>& aircraft)
        : ControlAllocatorBase<ControlSize, EffectorSize>(aircraft){
        std::cout << "DP_LP_ControlAllocator使用传入的aircraft对象给ControlAllocatorBase初始化aircraft成员后, 其余参数DP_LPCA_problem在这里初始化"<< std::endl;
        // 在此处用aircraft, generalizedMoment初始化 成员变量 DP_LPCA_problem 和 Pre_DP_LPCA_problem 
        // 线性规划数据
        DP_LPCA_problem.tol=1e-7;
        DP_LPCA_problem.itlim = 10;
        //填数据
        DP_LPCA_problem.inB[0]=0;
        DP_LPCA_problem.inB[1]=1;
        DP_LPCA_problem.inB[2]=3;

        DP_LPCA_problem.e[0] = true;
        DP_LPCA_problem.e[1] = true;
        DP_LPCA_problem.e[2] = false;
        DP_LPCA_problem.e[3] = true;
        DP_LPCA_problem.e[4] = true;
        
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            float temp=0;
            for(int j=0; j<DP_LPCA_problem.n-1; ++j)
            {
                DP_LPCA_problem.A[i][j] =aircraft.controlEffectMatrix[i][j];
                temp +=-aircraft.controlEffectMatrix[i][j]*aircraft.lowerLimits[j];
            }
            DP_LPCA_problem.A[i][DP_LPCA_problem.n-1] =0;
            DP_LPCA_problem.b[i] = temp;
        }
        for(int i=0; i<DP_LPCA_problem.n-1; ++i)
        {
            DP_LPCA_problem.c[i] =0;
        }
        DP_LPCA_problem.c[DP_LPCA_problem.n-1] =-1;
        for(int i=0; i<DP_LPCA_problem.n; ++i)
        {
            DP_LPCA_problem.h[i] =aircraft.upperLimits[i]-aircraft.lowerLimits[i];
        }
        DP_LPCA_problem.h[DP_LPCA_problem.n-1]=upper_lam;
        
        LPsolverForAC.problem=DP_LPCA_problem;
    }

    // 设置算法参数函数

    // 析构函数
    ~DP_LP_ControlAllocator() {
        // 如果有需要释放的资源，可以在这里添加代码
    }

    float* allocateControl(float input[ControlSize]) override {
        // 重写控制分配器函数
        // 实现控制分配算法
        float* EffectorCommand = new float[EffectorSize];
        // 算法实现
        // DP_LPCA（generalizedMoment, aircraft） 
        // DP_LPCA函数利用飞行器数据，将分配问题描述为DP_LP问题并用BoundedRevisedSimplex求解
        // 使用模版函数result = BoundedRevisedSimplex(problem); 
        //=======================
        bool flag=false;
        for(int i=0;i<ControlSize;++i){
            if(std::abs(input[i]) > DP_LPCA_problem.tol)
            {
                flag = true; // Check this condition
                break;
            }
        }
        if(!flag){
            for(int i=0;i<EffectorSize;++i){
                EffectorCommand[i]=0;
                // this->generalizedMoment[i] =0; //error do not add this
            }
            // std::cout << "return 0"<< std::endl;
            return EffectorCommand;
        }
        //=======================
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            DP_LPCA_problem.A[i][DP_LPCA_problem.n-1] =-input[i];
            this->generalizedMoment[i] = input[i]; // just record.
        }
        
        auto result = BoundedRevisedSimplex(DP_LPCA_problem);
        // 使用结果
        // result.y0, result.inB, result.e, result.errout
        int err = 0;
        float xout[DP_LPCA_problem.n];
        for(int i=0;i<DP_LPCA_problem.n;++i){
            xout[i]=0;
        }
        for(int i=0;i<ControlSize;++i){
            xout[result.inB[i]]=result.y0[i];
        }
        for(int i=0;i<DP_LPCA_problem.n;++i){
            if(!result.e[i]){
                xout[i]=-xout[i]+DP_LPCA_problem.h[i];
            }
        }
        if(result.iters>=DP_LPCA_problem.itlim){
            err = 3;
            std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
            for (int i = 0; i < ControlSize; ++i) {
                std::cout << this->generalizedMoment[i] << std::endl;
            }
        }
        if(result.errout)
        {
            err = 1;
            std::cout << "Solver error"<< std::endl;
            for (int i = 0; i < ControlSize; ++i) {
                std::cout << this->generalizedMoment[i] << std::endl;
            }
        }
        for(int i=0;i<EffectorSize;++i){
            EffectorCommand[i]=xout[i]+this->aircraft.lowerLimits[i];
        }
        if(xout[EffectorSize]>1){
            for(int i=0;i<EffectorSize;++i){
                EffectorCommand[i]/=xout[EffectorSize];
            }
        }
        return EffectorCommand;
    }
    float* allocateControl_bases_solver(float input[ControlSize]){
        // 重写控制分配器函数
        // 实现控制分配算法
        float* EffectorCommand = new float[EffectorSize];
        // 算法实现
        // DP_LPCA（generalizedMoment, aircraft） 
        // DP_LPCA函数利用飞行器数据，将分配问题描述为DP_LP问题并用BoundedRevisedSimplex求解
        // 使用模版函数result = BoundedRevisedSimplex(problem);  
        for (int i = 0; i < ControlSize; ++i) {
            LPsolverForAC.problem.A[i][EffectorSize]=-input[i];
            this->generalizedMoment[i] = input[i];
        }
        auto result = LPsolverForAC.solve();
        // 使用结果
        // result.y0, result.inB, result.e, result.errout
        int err = 0;
        float xout[LPsolverForAC.problem.n];
        for(int i=0;i<LPsolverForAC.problem.n;++i){
            xout[i]=0;
        }
        for(int i=0;i<ControlSize;++i){
            xout[result.inB[i]]=result.y0[i];
        }
        for(int i=0;i<LPsolverForAC.problem.n;++i){
            if(!result.e[i]){
                xout[i]=-xout[i]+LPsolverForAC.problem.h[i];
            }
        }
        if(result.iters>=LPsolverForAC.problem.itlim){
            err = 3;
            std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
        }
        if(result.errout)
        {
            err = 1;
            std::cout << "Solver error"<< std::endl;
        }
        for(int i=0;i<EffectorSize;++i){
            EffectorCommand[i]=xout[i]+this->aircraft.lowerLimits[i];
            // EffectorCommand[i]=0;
        }
        if(xout[EffectorSize]>1){
            for(int i=0;i<EffectorSize;++i){
                EffectorCommand[i]/=xout[EffectorSize];
            }
        }
        return EffectorCommand;
    }
    // 其他成员函数和成员变量定义
    float generalizedMoment[ControlSize]; // 构造函数设置
    // 线性规划相关
    LinearProgrammingProblem<ControlSize, EffectorSize+1> DP_LPCA_problem;// 提前设置 using aircraft
    LinearProgrammingProblem<ControlSize, EffectorSize+4> Pre_DP_LPCA_problem;// 提前设置 using aircraft
    LinearProgrammingSolverForAC<ControlSize, EffectorSize+1> LPsolverForAC;
    float upper_lam = 1e4;
    
};

// template <int ControlSize, int EffectSize>
// class XX_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectSize> {}