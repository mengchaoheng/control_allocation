
#include <matrix/math.hpp>
#include <iostream>

using namespace matrix;
// Add the min_user function definition here
template<typename Type, size_t M, size_t N>
inline void min_user(const Matrix<Type, M, N> &x, Type &x_min, size_t &x_index)
{
    const Matrix<Type, M, N> &self = x;
    x_min = self(0, 0);
    x_index = 0;

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (self(i, j) < x_min) {
                x_min = self(i, j);
                x_index = i; // Assuming that you want the index along the row
            }
        }
    }
}
// 计算 rho 的函数, for DPscaled_LPCA
const int SIZE_ydt = 3; // 假设 ydt 是一个包含 5 个元素的一维数组
const int SIZE_Bt_row = 3; // 假设 Bt 是一个 5x5 的二维数组
const int SIZE_Bt_col = 4;
inline float calculateRho(float ydt[], float u[], float Bt[][SIZE_Bt_col], float tol) {
    float numerator = 0.0f;
    float denominator = 0.0f;
    float ydt_T_Bt[SIZE_Bt_col];
    for (int j = 0; j < SIZE_Bt_col; ++j) {
        ydt_T_Bt[j] = 0;
        for (int k = 0; k < SIZE_ydt; ++k) { //or < SIZE_Bt_row
            ydt_T_Bt[j] += ydt[k] * Bt[k][j];
        }
    }
    for (int k = 0; k < SIZE_Bt_col; ++k) {
        numerator += ydt_T_Bt[k] * u[k];
    }
    // 计算 ydt 的模的平方
    for (int i = 0; i < SIZE_ydt; ++i) {
        denominator += ydt[i] * ydt[i];
    }
    // 避免除以零  // ydt的模不会很小
    if (std::abs(denominator) < tol ) {
        // std::cerr << "Error: Division by zero." << std::endl;
        return 0.0f;
    }
    // 计算 rho
    return numerator / denominator;
}

// 函数用于计算两个正整数集合的差
inline void setdiff(int setA[], int sizeA, int setB[], int sizeB, int result[]) {
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
inline int* generateSequence(int i, int n) {
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
    // 赋值运算符
    LinearProgrammingProblem& operator=(const LinearProgrammingProblem& other) {
        if (this != &other) {
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
        return *this;
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


// 定义函数模板，修正单纯形算法实现。要求A行满秩，求解问题前已知inB和e，即需要找到一个初始基本可行解开始算法迭代。e=0表示该初始解在上限h上，否则就是0
// see A.6.3 Simplex Method[1]. this function reimplement the simplxuprevsol of this book:
// [1] W. Durham, K. A. Bordignon, and R. Beck, Aircraft control allocation. none: John Wiley & Sons, 2017.
// and you can download the code on: https://github.com/mengchaoheng/control_allocation.git
template<int M, int N>
LinearProgrammingResult<M, N> BoundedRevisedSimplex(LinearProgrammingProblem<M, N> problem) {
    // Bounded Revised Simplex

    // function [yout, inBout,eout, itout,errout] = simplxuprevsol(A,ct,b,inB,inD,h,e,m,n,itlim)

    // Solves the linear program:
    //         minimize c'y 
    //         subject to 
    //         Ay = b
    //         0<= y <= h

    // Inputs: 
    //         A [m,n]   = lhs Matrix of equaltity constraints
    //         ct [1,n]  = transpose of cost vector
    //         b [m,1]   = rhs vector for equality constraint
    //         inB [m]   = Vector of indices of unknowns in the initial basic set
    //         inD [n-m] = Vector of indices of unknowns not in the initial basic set
    //         h[n,1]    = Upper Bound for unknowns
    //         e[n,1]    = Sign for unknown variables (+ lower bound, - upper bound)
    // Optional inputs:
    //         m,n       = number of constraints, unknowns (Opposite standard
    //                     CA convention
    //         itlim     = Upper bound on the allowed iterations

    // Outputs:
    //         yout[n,1]  = Optimal output variable
    //         inBout     = indices of Basic vectors in output
    //         eout       = sign associate with output unknowns
    //         itout      = number of iterations remaining out of itlim
    //         errout     = Flag (=true) if unbounded is set

    // Modification History
    // 2002      Roger Beck  Original
    // 8/2014    Roger Beck  Update for use
    // 9/2014    Roger Beck  Added anti-cycling rule
    // 4/2024    Meng ChaoHeng
    // std::cout << "BoundedRevisedSimplex loop"<< std::endl; 
    // std::cout << "M"<<M<< std::endl; 
    // std::cout << "N"<<N<< std::endl; 

    LinearProgrammingResult<M, N> result;
    // 使用 problem.inB, problem.inD, problem.itlim, problem.A, problem.b, problem.c, problem.h, problem.e, problem.tol
    const int n_m=N-M;
    // Index list for non-basic variables, that is 1 2 3 4 ... n
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
    // Index list for all variables 
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
    // Partition A, we have inD, which the element in ind_all but not in inB
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

    
    //djust signs problem if variables are initialized at upper bounds.
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
    // std::cout << "e: [";
    // for (size_t i = 0; i < N; ++i) {
    //     std::cout << problem.e[i];
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
    //===============use inB and inD to inital==============
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<M; ++j)
        {
            A_inB(i,j)=problem.A[i][problem.inB[j]];
        }
        c_inB(i)=problem.c[problem.inB[i]];
    }
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j<N-M; ++j)
        {
            A_inD(i,j)=problem.A[i][problem.inD[j]];
        }
    }
    for(int i=0; i<N-M; ++i)
    {
        c_inD(i)=problem.c[problem.inD[i]];
    }
    matrix::Vector<float, M> b_vec(problem.b);
    
    // std::cout << "A_inB:";
    // A_inB.print();
    // std::cout << "A_inD:";
    // A_inD.print();
    // std::cout << "c_inB:";
    // c_inB.print();
    // std::cout << "c_inD:";
    // c_inD.print();
    // std::cout << "b_vec:";
    // b_vec.print();
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
    //==============================

    //  %Initial Solution
    matrix::LeastSquaresSolver<float, M,M> LSsolver0(A_inB);
    matrix::Vector<float, M> y0 = LSsolver0.solve(b_vec);
    // std::cout << "y0:";
    // y0.print();
    // Initialize Loop Termination Conditionss
    bool done = false;
    bool unbounded = false;
    int iters =0;
     while ((!done  || !unbounded ) && (iters <= problem.itlim))
    {
        iters = iters+1;
        // Calculate transpose of relative cost vector based on current basis
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
        // Find minimum relative cost
        min_user(rdt.transpose(), minr, qind);
        // std::cout << "minr:"<<minr<<std::endl;
        // std::cout << "qind:"<<qind<<std::endl;
        if(minr >=0)  // If all relative costs are positive then the solution is optimal. have to compare with 0 !
        { 
            done = true;
            break;
        }
        int qel = problem.inD[qind];  // Unknown to Enter the basis minimizes relative cost
        // A_qel(0)=problem.A[0][qel];
        // A_qel(1)=problem.A[1][qel];
        // A_qel(2)=problem.A[2][qel];
        for(int i=0;i<problem.m;++i){
            A_qel(i)=problem.A[i][qel];
        }
        // std::cout << "qel:"<<qel<<std::endl;
        // yq=inv(A_inB)* A_qel; 
        matrix::LeastSquaresSolver<float, M,M> LSsolver1(A_inB);
        yq = LSsolver1.solve(A_qel); // Vector to enter in terms of the current Basis vector
        // std::cout << "yq:";
        // yq.print();
         // Check wether all the abs of yq[i] is greater than tol.
        bool flag=false;
        for(int i=0;i<M;++i){
            if(std::abs(yq(i)) > problem.tol) 
            {
                flag = true; 
                break;
            }
        }
        if(!flag)
        {
            unbounded = true; // Check this condition
            std::cout << "simplex loop Solution is unbounded"<< std::endl; 
            break;
        }
        // Compute ratio how much each current basic variable will have to move for the entering variable.
        // careful here
        float hinB[M];
        for(int i=0;i<M;++i)
        {
            if(std::abs(yq(i))>problem.tol)
            {
                rat(i)=y0(i)/yq(i);
                hinB[i]=problem.h[problem.inB[i]];
                // If yq < 0 then increasing variable when it leaves the basis will minimize cost
                if(yq(i)<0 ) // have to be compare with 0!!!
                {
                    // std::cout << "yq(i)<0 i:"<< i<< std::endl;                  
                    rat(i)-=hinB[i]/yq(i);
                }
            }
            else // If an element yq ~=0 then it doesn't change for the entering variable and shouldn't be chosen
            {
                rat(i)=INFINITY;
            }
        }
        // std::cout << "rat:";
        // rat.print();
         // Variable to exit is moving to its minimum value--Note that min_user returns the lowest index minimum
        float minrat=rat(0);
        size_t p=0;
        min_user(rat, minrat, p);
        // std::cout << "minrat:"<<minrat<<std::endl;
        // std::cout << "p:"<<p<<std::endl;

        // If the minimum ratio is zero, then the solution is degenerate and the entering
        // variable will not change the basis---invoke Bland's selection rule to avoid
        // cycling.
        if (std::abs(minrat) <= problem.tol)
        {
            // Find negative relative cost
            // std::cout << " Find negative relative cost "<< std::endl; 
            for(int i=0;i<N-M;++i)
            {
                // std::cout << "rdt(0,i):"<<rdt(0,i)<<std::endl; 
                // indm is the index of rdt < 0, qind is the fisrt one.
                if(rdt(0,i)<0){ // Note that since minr <0 indm is not empty 
                    qind=nind[i];
                    qel = problem.inD[qind];// Unknown to Enter the basis is first indexed to avoid cycling
                    // std::cout << "qind:"<<qind<<std::endl;
                    // std::cout << "qel:"<<qel<<std::endl;
                    break;
                }
            }
            // A_qel(0)=problem.A[0][qel];
            // A_qel(1)=problem.A[1][qel];
            // A_qel(2)=problem.A[2][qel];
            for(int i=0;i<problem.m;++i){
                A_qel(i)=problem.A[i][qel];
            }
            // yq=inv(A_inB)* A_qel;
            matrix::LeastSquaresSolver<float, M,M> LSsolver2(A_inB);
            yq = LSsolver2.solve(A_qel); // Vector to enter in terms of the current Basis vector
            // std::cout << "yq:";
            // yq.print();
            bool flag1=false;
            for(int i=0;i<M;++i){
                if(std::abs(yq(i)) > problem.tol)
                {
                    flag1 = true; 
                    break;
                }
            }
            if(!flag1)
            {
                unbounded = true; // Check this condition
                std::cout << "simplex loop Solution is unbounded"<< std::endl; 
                break;
            }
            // Recompute rations and determine variable to leave
            float hinB1[M];
            for(int i=0;i<M;++i)
            {
                hinB1[i]=problem.h[problem.inB[i]];
                if(std::abs(yq(i))>problem.tol) 
                {
                    rat(i)=y0(i)/yq(i);
                    if(yq(i)<0) // If yq < 0 then increasing variable when it leaves the basis will minimize cost
                    {                    
                        rat(i)-=hinB1[i]/yq(i);
                    }
                }
                else
                { 
                    rat(i)=INFINITY; // If an element yq ~=0 then it doesn't change for the entering variable and shouldn't be chosen
                }
            }
            // std::cout << "rat:";
            // rat.print();
            // Variable to exit is moving to its minimum value--Note that min_user returns the lowest index minimum
            minrat=rat(0);
            p=0;
            min_user(rat, minrat, p);
            // std::cout << "minrat:"<<minrat<<std::endl;
            // std::cout << "p:"<<p<<std::endl;
        }
        // Maintain the bounded simplex as only having lower bounds by recasting any variable that needs to move to its opposite bound.
        // std::cout << "problem.h[qel]:"<<problem.h[qel]<<std::endl;
        // std::cout << "minrat-problem.h[qel]:"<<minrat-problem.h[qel]<<std::endl;
        if (minrat >= problem.h[qel])
        {
            // Case 1: Entering variable goes to opposite bound and current basis is maintained
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
            // Case 2: Leaving variable returns to lower bound (0)	
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
            // Case 2: Leaving variable moves to upper bound	
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
        //  Compute new Basic solution;
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
        // std::cout << "e: [";
        // for (size_t i = 0; i < N; ++i) {
        //     std::cout << problem.e[i];
        //     if (i < N- 1) {
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
    Aircraft(const float (&controlEffectMatrixInit)[ControlSize][EffectorSize], 
             const float (&upperLimitsInit)[EffectorSize], 
             const float (&lowerLimitsInit)[EffectorSize]) {
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
    Aircraft(const float (&upperLimitsInit)[EffectorSize], 
             const float (&lowerLimitsInit)[EffectorSize]) {
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
    Aircraft(const float (&controlEffectMatrixInit)[ControlSize][EffectorSize], const float& lowerBound, const float& upperBound) : lower(lowerBound), upper(upperBound){
        // 使用传入的初始化参数和飞行器
        // std::cout << "使用传入的controlEffectMatrixInit, lower, upper初始化参数和飞行器"<< std::endl; 
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
    float lower;
    float upper;

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
    ControlAllocatorBase(const Aircraft<ControlSize, EffectorSize>& ac) 
        : aircraft(ac) { 
        // 使用传入的aircraft对象初始化aircraft成员
        // std::cout << "ControlAllocatorBase使用传入的aircraft对象初始化aircraft成员"<< std::endl;
    }

    virtual void allocateControl(float input[ControlSize], float output[EffectorSize], int& err) = 0;

    // 其他数学函数和成员变量定义
    Aircraft<ControlSize, EffectorSize> aircraft; // 构造函数设置
    
};
// 控制分配类模板
template <int ControlSize, int EffectorSize>
class DP_LP_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectorSize> {
private:
    // 添加算法设置参数
public:
    // 构造函数, 利用aircraft 预设置LinearProgrammingProblem
    // 构造函数
    DP_LP_ControlAllocator(const Aircraft<ControlSize, EffectorSize>& ac)
        : ControlAllocatorBase<ControlSize, EffectorSize>(ac){
        // std::cout << "DP_LP_ControlAllocator使用传入的aircraft对象给ControlAllocatorBase初始化aircraft成员后, 其余参数DP_LPCA_problem在这里初始化"<< std::endl;
        // 在此处用aircraft, generalizedMoment初始化 成员变量 DP_LPCA_problem 和 Pre_DP_LPCA_problem 
        // 线性规划数据
        //=====================================DP_LPCA_problem================================
        DP_LPCA_problem.tol=1e-7;
        DP_LPCA_problem.itlim = 10;
        //填数据
        
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            float temp=0;
            for(int j=0; j<DP_LPCA_problem.n-1; ++j)
            {
                DP_LPCA_problem.A[i][j] = this->aircraft.controlEffectMatrix[i][j];
                temp += -this->aircraft.controlEffectMatrix[i][j]*this->aircraft.lowerLimits[j];
            }
            DP_LPCA_problem.A[i][DP_LPCA_problem.n-1] = 0;
            DP_LPCA_problem.b[i] = temp;
        }
        for(int i=0; i<DP_LPCA_problem.n-1; ++i)
        {
            DP_LPCA_problem.c[i] = 0;
        }
        DP_LPCA_problem.c[DP_LPCA_problem.n-1] = -1;
        for(int i=0; i<DP_LPCA_problem.n-1; ++i)
        {
            DP_LPCA_problem.h[i] = this->aircraft.upperLimits[i]-this->aircraft.lowerLimits[i];
        }
        DP_LPCA_problem.h[DP_LPCA_problem.n-1] = upper_lam;
        
        //==================================PreDP_LPCA_problem================================
        Pre_DP_LPCA_problem.tol=1e-7;
        Pre_DP_LPCA_problem.itlim = 10;
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            for(int j=0; j<DP_LPCA_problem.n; ++j)
            {
                Pre_DP_LPCA_problem.A[i][j] = DP_LPCA_problem.A[i][j];
            }
            Pre_DP_LPCA_problem.b[i] = DP_LPCA_problem.b[i]; // the same as DP_LPCA_problem

        }
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            Pre_DP_LPCA_problem.A[i][i + DP_LPCA_problem.n] = (DP_LPCA_problem.b[i] > 0) ? 1 : -1; // sb = 2*(b > 0)-1; Ai = [A diag(sb)]; 
        }
        for(int i=0; i<DP_LPCA_problem.n; ++i)
        {
            Pre_DP_LPCA_problem.c[i] =0;
        }
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            Pre_DP_LPCA_problem.c[i+DP_LPCA_problem.n] = 1;
        }
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            Pre_DP_LPCA_problem.inB[i] = DP_LPCA_problem.n+i;
        }
        for(int i=0; i<DP_LPCA_problem.m+DP_LPCA_problem.n; ++i)
        {
            Pre_DP_LPCA_problem.e[i] = true;
        }
        for(int i=0; i<DP_LPCA_problem.n; ++i)
        {
            Pre_DP_LPCA_problem.h[i] = DP_LPCA_problem.h[i];
        }
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            Pre_DP_LPCA_problem.h[i+DP_LPCA_problem.n] = 2*std::abs(DP_LPCA_problem.b[i]);
        }
        //================================== DPscaled_LPCA_problem ================================
        // std::cout << "DPscaled_LPCA_problem"<< std::endl; 
        DPscaled_LPCA_problem.tol=1e-7;
        DPscaled_LPCA_problem.itlim = 10;
        float yd[3]={0.1,0.2,-0.1}; // random value for inital.
        float my=yd[0]; // 存储最大的绝对值
        int iy=0; // 存储最大绝对值的索引
        for (int i = 0; i < ControlSize; ++i) {
            float absValue = std::abs(yd[i]); // 计算 yd 中第 i 个元素的绝对值
            if (absValue > my) { // 如果当前绝对值大于 my，则更新 my 和 iy
                my = absValue;
                iy = i;
            }
        }
        // std::cout << "DPscaled_LPCA_problem iy: "<< iy <<std::endl; 
        // std::cout << "DPscaled_LPCA_problem my: "<< my <<std::endl; 
        // copy firstly !!!
        float Bt[ControlSize][EffectorSize];
        float ydt[ControlSize];
        
        for(int i=0;i<ControlSize;++i){
            for (int j = 0; j <EffectorSize; ++j) {
                Bt[i][j] = this->aircraft.controlEffectMatrix[i][j];
            }
        }
        for(int i=0;i<ControlSize;++i){
            ydt[i]=yd[i];
        }

        // move
        for(int j=0;j<EffectorSize;++j){
            Bt[0][j] = this->aircraft.controlEffectMatrix[iy][j];
            for (int i = iy; i >0; i--) {
                Bt[i][j] = this->aircraft.controlEffectMatrix[i-1][j];
            }
        }
        ydt[0] = yd[iy];
        for (int j = iy; j >0; j--) {
            ydt[j] = yd[j-1];
        }
        // swap 2 3 row of ydt
        float ydt2=ydt[1];
        ydt[1]=ydt[2];
        ydt[2]=ydt2;
        // swap 2 3 col of Bt
        for(int i=0;i<EffectorSize;++i){
            float Bt2=Bt[1][i];
            Bt[1][i]=Bt[2][i];
            Bt[2][i]=Bt2;
        }
        
        // std::cout << "aircraft.controlEffectMatrix:" << std::endl;
        // for (size_t i = 0; i < ControlSize; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << this->aircraft.controlEffectMatrix[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "Bt:" << std::endl;
        // for (size_t i = 0; i < ControlSize; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << Bt[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // M = [ydt(2:ControlSize) -ydt(1)*eye(ControlSize-1)];
        float M[ControlSize-1][ControlSize];
        // M[0][0]=ydt[1];
        // M[1][0]=ydt[2];
        // M[0][1]=-ydt[0];
        // M[1][1]=0;
        // M[0][2]=0;
        // M[1][2]=-ydt[0];
        // or 
        // 将 M 的所有元素初始化为0，并同时填充 M 的第一列和对角线元素
        for (int i = 0; i < ControlSize - 1; ++i) {
            for (int j = 0; j < ControlSize; ++j) {
                if (j == 0) {
                    // 填充 M 的第一列
                    M[i][0] = ydt[i + 1];
                } else if (j == i + 1) {
                    // 填充对角线元素
                    M[i][j] = -ydt[0];
                } else {
                    // 初始化其他元素为0
                    M[i][j] = 0.0f;
                }
            }
        }

        // std::cout << "M:" << std::endl;
        // for (size_t i = 0; i < ControlSize-1; ++i) {
        //     for (size_t j = 0; j < ControlSize; ++j) {
        //         std::cout << M[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        for (int i = 0; i < ControlSize-1; ++i) {
            for (int j = 0; j < EffectorSize; ++j) {
                DPscaled_LPCA_problem.A[i][j] = 0;
                for (int k = 0; k < ControlSize; ++k) {
                    DPscaled_LPCA_problem.A[i][j] += M[i][k] * Bt[k][j];
                }
            }
        }
        // std::cout << "DPscaled_LPCA_problem.A:" << std::endl;
        // for (size_t i = 0; i < ControlSize-1; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << DPscaled_LPCA_problem.A[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        for(int i=0; i<ControlSize-1; ++i)
        {
            float temp=0;
            for(int j=0; j<EffectorSize; ++j)
            {
                temp += -DPscaled_LPCA_problem.A[i][j]*this->aircraft.lowerLimits[j];
            }
            DPscaled_LPCA_problem.b[i] = temp;
        }
        for (int i = 0; i < EffectorSize; ++i) {
            DPscaled_LPCA_problem.c[i] = 0;
            for (int j = 0; j < ControlSize; ++j) {
                DPscaled_LPCA_problem.c[i] += -Bt[j][i] * ydt[j];
            }
        }
        for (int i = 0; i < EffectorSize; ++i) {
            DPscaled_LPCA_problem.h[i] = this->aircraft.upperLimits[i]-this->aircraft.lowerLimits[i];
        }
        // std::cout << "DPscaled_LPCA_problem b: [";
        // for (size_t i = 0; i < ControlSize-1; ++i) {
        //     std::cout << DPscaled_LPCA_problem.b[i];
        //     if (i < ControlSize-1 - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;


        // std::cout << "DPscaled_LPCA_problem c: [";
        // for (size_t i = 0; i <EffectorSize; ++i) {
        //     std::cout << DPscaled_LPCA_problem.c[i];
        //     if (i <EffectorSize - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "DPscaled_LPCA_problem h: [";
        // for (size_t i = 0; i < EffectorSize; ++i) {
        //     std::cout << DPscaled_LPCA_problem.h[i];
        //     if (i < EffectorSize- 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        //==================================Pre_DPscaled_LPCA_problem================================
        // std::cout << "Pre_DPscaled_LPCA_problem"<< std::endl; 
        Pre_DPscaled_LPCA_problem.tol=1e-7;
        Pre_DPscaled_LPCA_problem.itlim = 10;
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            for(int j=0; j<DPscaled_LPCA_problem.n; ++j)
            {
                Pre_DPscaled_LPCA_problem.A[i][j] = DPscaled_LPCA_problem.A[i][j];
            }
            Pre_DPscaled_LPCA_problem.b[i] = DPscaled_LPCA_problem.b[i]; // the same as DPscaled_LPCA_problem

        }
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            Pre_DPscaled_LPCA_problem.A[i][i + DPscaled_LPCA_problem.n] = (DPscaled_LPCA_problem.b[i] > 0) ? 1 : -1; // sb = 2*(b > 0)-1; Ai = [A diag(sb)]; 
        }
        for(int i=0; i<DPscaled_LPCA_problem.n; ++i)
        {
            Pre_DPscaled_LPCA_problem.c[i] =0;
        }
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            Pre_DPscaled_LPCA_problem.c[i+DPscaled_LPCA_problem.n] = 1;
        }
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            Pre_DPscaled_LPCA_problem.inB[i] = DPscaled_LPCA_problem.n+i;
        }
        for(int i=0; i<DPscaled_LPCA_problem.m+DPscaled_LPCA_problem.n; ++i)
        {
            Pre_DPscaled_LPCA_problem.e[i] = true;
        }
        for(int i=0; i<DPscaled_LPCA_problem.n; ++i)
        {
            Pre_DPscaled_LPCA_problem.h[i] = DPscaled_LPCA_problem.h[i];
        }
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            Pre_DPscaled_LPCA_problem.h[i+DPscaled_LPCA_problem.n] = 2*std::abs(DPscaled_LPCA_problem.b[i]);
        }

        // std::cout << "Pre_DPscaled_LPCA_problem A:" << std::endl;
        // for (int i = 0; i < Pre_DPscaled_LPCA_problem.m; ++i) {
        //     for (int j = 0; j < Pre_DPscaled_LPCA_problem.n; ++j) {
        //         std::cout << Pre_DPscaled_LPCA_problem.A[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "Pre_DPscaled_LPCA_problem b: [";
        // for (int i = 0; i < Pre_DPscaled_LPCA_problem.m; ++i) {
        //     std::cout << Pre_DPscaled_LPCA_problem.b[i];
        //     if (i < Pre_DPscaled_LPCA_problem.m - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        // std::cout << "Pre_DPscaled_LPCA_problem c: [";
        // for (int i = 0; i < Pre_DPscaled_LPCA_problem.n; ++i) {
        //     std::cout << Pre_DPscaled_LPCA_problem.c[i];
        //     if (i < Pre_DPscaled_LPCA_problem.n - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        // std::cout << "Pre_DPscaled_LPCA_problem h: [";
        // for (int i = 0; i < Pre_DPscaled_LPCA_problem.n; ++i) {
        //     std::cout << Pre_DPscaled_LPCA_problem.h[i];
        //     if (i < Pre_DPscaled_LPCA_problem.n- 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
    }

    // 设置算法参数函数

    // 析构函数
    ~DP_LP_ControlAllocator() {
        // 如果有需要释放的资源，可以在这里添加代码
    }

    // 最初测试算法使用本函数，由于DP_LPCA和4片舵涵道的特殊性，初始基本解是可以预先确定且不变的。所以可以省略第一步寻找基本初始解。后续可以将本函数改为切换使用DP_LPCA和DPscaled_LPCA。
    void allocateControl(float input[ControlSize], float output[EffectorSize], int& err) override {
        // 重写控制分配器函数
        // 实现控制分配算法
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
                output[i]=0;
            }
            // std::cout << "return 0"<< std::endl;
            return;
        }
        //=======================
        //===========just for df4, we alway have to calc this by a new problem================
        // we can call BoundedRevisedSimplex direct in allocationControl
        DP_LPCA_problem.inB[0]=0;
        DP_LPCA_problem.inB[1]=1;
        DP_LPCA_problem.inB[2]=3;
        DP_LPCA_problem.e[0] = true;
        DP_LPCA_problem.e[1] = true;
        DP_LPCA_problem.e[2] = false;
        DP_LPCA_problem.e[3] = true;
        DP_LPCA_problem.e[4] = true;
        //=======================
        // update A b h
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            float temp=0;
            for(int j=0; j<DP_LPCA_problem.n-1; ++j)
            {
                DP_LPCA_problem.A[i][j] = this->aircraft.controlEffectMatrix[i][j];
                temp += -this->aircraft.controlEffectMatrix[i][j]*this->aircraft.lowerLimits[j];
            }
            DP_LPCA_problem.A[i][DP_LPCA_problem.n-1] = -input[i];
            this->generalizedMoment[i] = input[i]; // just record.
            DP_LPCA_problem.b[i] = temp;
        }
        for(int i=0; i<DP_LPCA_problem.n-1; ++i)
        {
            DP_LPCA_problem.h[i] = this->aircraft.upperLimits[i]-this->aircraft.lowerLimits[i];
        }
        
        auto result = BoundedRevisedSimplex(DP_LPCA_problem);
        // 使用结果
        // result.y0, result.inB, result.e, result.errout
        // int err = 0;
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
            // std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
            // for (int i = 0; i < ControlSize; ++i) {
            //     std::cout << this->generalizedMoment[i] << std::endl;
            // }
        }
        if(result.errout)
        {
            err = 1;
            // std::cout << "Solver error"<< std::endl;
            // for (int i = 0; i < ControlSize; ++i) {
            //     std::cout << this->generalizedMoment[i] << std::endl;
            // }
        }
        for(int i=0;i<EffectorSize;++i){
            output[i]=xout[i]+this->aircraft.lowerLimits[i];
        }
        if(xout[EffectorSize]>1){
            for(int i=0;i<EffectorSize;++i){
                output[i]/=xout[EffectorSize];
            }
        }
        return;
    }
    
    // To find an initial condition, many linear programming solvers treat the solution in two phases. Phase one solves a specially constructed problem designed to yield a basic feasible solution that is used to initialize the original problem in phase two.
    // So we have DP_LPCA and DPscaled_LPCA
    void DP_LPCA(float input[ControlSize], float output[EffectorSize], int& err){
        // Direction Preserving Control Allocation Linear Program

        // function [u, errout] = DP_LPCA(yd,B,uMin,uMax,itlim,upper_lam);

        // Solves the control allocation problem while preserving the
        // objective direction for unattainable commands. The solution
        // is found by solving the problem,
        // min -lambda,
        // s.t. B*u = lambda*yd, uMin<=u<=uMax, 0<= lambda <=1

        // For yd outside the AMS, the solution returned is that the
        // maximum in the direction of yd.

        // For yd strictly inside the AMS, the solution achieves
        // Bu=yd and m-n controls will be at their limits; but there
        // is no explicit preference to which solution will be 
        // returned. This limits the usefulness of this routine as
        // a practical allocator unless preferences for attainable solutions
        // are handled externally.

        // (For derivation of a similar formulation see A.1.2 and A.2.3 in the
        // text)


        // Inputs:
        //         yd [n]    = Desired objective
        //         B [n,m]   = Control Effectiveness matrix
        //         uMin[m,1] = Lower bound for controls
        //         uMax[m,1] = Upper bound for controls
        //         itlim     = Number of allowed iterations limit
        //                         (Sum of iterations in both branches)

        // Outputs:
        //         u[m,1]     = Control Solution
        //         errout     = Error Status code
        //                         0 = found solution
        //                         <0 = Error in finding initial basic feasible solution
        //                         >0 = Error in finding final solution
        //                         -1,1 = Solver error (unbounded solution)
        //                         -2   = Initial feasible solution not found
        //                         -3,3 = Iteration limit exceeded
        //         itlim      = Number of iterations remaining after solution found
        //         upper_lam  = the upper limit of lambda

        // Calls:
        //         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)

        // Notes:
        // If errout ~0 there was a problem in the solution. %

        // Error code < 0 implies an error in the initialization and there is no guarantee on
        // the quality of the output solution other than the control limits.
        // Error code > 0 for errors in final solution--B*u is in the correct direction and has
        // magnitude < yd, but B*u may not equal yd (for yd attainable)
        // or be maximized (for yd unattainable)

        // Modification History
        // 2002      Roger Beck  Original (DPcaLP8.m)
        // 8/2014    Roger Beck  Update for use in text
        // 4/2024    Meng ChaoHeng  Implement in cpp 
        // DP_LPCA函数利用飞行器数据，将分配问题描述为DP_LP问题并用BoundedRevisedSimplex求解
        //=======================
        // std::cout << "DP_LPCA"<< std::endl; 
        // Initialize error code to zero
        // int err = 0;
        // Figure out how big the problem is (use standard CA definitions for m & n)
        // but in here we use [m,k] = size(B) instead of [n,m] = size(B) in matlab. just for adapt to BoundedRevisedSimplex
        // we use [m,n] = size(A) in BoundedRevisedSimplex, that is, k + 1 = n.
        // Check to see if yd == 0
        // May want to adjust the tolerance to improve numerics of later steps
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
                output[i]=0;
            }
            // std::cout << "return 0"<< std::endl;
            return;
        }
        //=======================  
        // We inital this problem in constructor.
        // Construct an LP using scaling parameter to enforce direction preserving
        // To find Feasible solution construct problem with appended slack variables
        // ref. is A.6.4 Initialization of the Simplex Algorithm of <Aircraft control allocation>

        // now we update the problem by input data.
        // update A b h
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            float temp=0;
            for(int j=0; j<DP_LPCA_problem.n-1; ++j)
            {
                DP_LPCA_problem.A[i][j] = this->aircraft.controlEffectMatrix[i][j];
                temp += -this->aircraft.controlEffectMatrix[i][j]*this->aircraft.lowerLimits[j];
            }
            DP_LPCA_problem.A[i][DP_LPCA_problem.n-1] = -input[i];
            this->generalizedMoment[i] = input[i]; // just record.
            DP_LPCA_problem.b[i] = temp;
        }
        for(int i=0; i<DP_LPCA_problem.n-1; ++i)
        {
            DP_LPCA_problem.h[i] = this->aircraft.upperLimits[i]-this->aircraft.lowerLimits[i];
        }
        // update sb(since b is update) Ai hi
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            for(int j=0; j<DP_LPCA_problem.n; ++j)
            {
                Pre_DP_LPCA_problem.A[i][j] = DP_LPCA_problem.A[i][j];
            }
            Pre_DP_LPCA_problem.b[i] = DP_LPCA_problem.b[i]; // the same as DP_LPCA_problem

        }
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            Pre_DP_LPCA_problem.A[i][i + DP_LPCA_problem.n] = (DP_LPCA_problem.b[i] > 0) ? 1 : -1; // sb = 2*(b > 0)-1; Ai = [A diag(sb)]; 
        }
        for(int i=0; i<DP_LPCA_problem.n; ++i)
        {
            Pre_DP_LPCA_problem.h[i] = DP_LPCA_problem.h[i];
        }
        for(int i=0; i<DP_LPCA_problem.m; ++i)
        {
            Pre_DP_LPCA_problem.h[i+DP_LPCA_problem.n] = 2*std::abs(DP_LPCA_problem.b[i]);
        }

        // Use Bounded Revised Simplex to find initial basic feasible point of original program
        auto result_init = BoundedRevisedSimplex(Pre_DP_LPCA_problem);

        // Check that Feasible Solution was found
        if(result_init.iters>=Pre_DP_LPCA_problem.itlim){
            err = -3;
            // std::cout << "Pre Too Many Iterations Finding Final Solution"<< std::endl; 
            // for (int i = 0; i < ControlSize; ++i) {
            //     std::cout << this->generalizedMoment[i] << std::endl;
            // }
        }
        for(int i=0;i<ControlSize;++i){
            if(result_init.inB[i]> EffectorSize) // DP_LPCA_problem is origin problem, k=DP_LPCA_problem.n-1 = EffectorSize 
            {
                // which mean inital basic index is out of the origin problem.
                err = -2;
                // std::cout << "Pre No Initial Feasible Solution found"<< std::endl; 
                // for (int k = 0; k < ControlSize; ++k) {
                //     std::cout << this->generalizedMoment[k] << std::endl;
                // }
                break;
            }
        }
        if(result_init.errout){
            err = -1;
            // std::cout << "Pre Solver error"<< std::endl; 
            // for (int i = 0; i < ControlSize; ++i) {
            //     std::cout << this->generalizedMoment[i] << std::endl;
            // }
        }
        // solve Pre_DP_LPCA_problem but proccess DP_LPCA_problem
        float xout[DP_LPCA_problem.n];
        for(int i=0;i<DP_LPCA_problem.n;++i){
            xout[i]=0;
        }
        if(err!=0) // Construct an incorrect solution to accompany error flags
        {
            // use result_init data
            // matlab: indv = inB1<=(k+1); xout(inB1(indv)) = y1(indv); % in matlab the index from 1 to k, but cpp is 0 to k-1
            for(int i=0;i<ControlSize;++i){
                if(result_init.inB[i] <= EffectorSize) 
                {
                    xout[result_init.inB[i]]=result_init.y0[i];
                }
            }
            // xout(~e1(1:k+1)) = -xout(~e1(1:k+1))+h(~e1(1:k+1));
            for(int i=0;i<DP_LPCA_problem.n;++i){
                if(!result_init.e[i]){
                    xout[i] = -xout[i] + DP_LPCA_problem.h[i];
                }
            }
        }
        else //No Error continue to solve problem
        { 
            // Solve using initial problem from above
            // Construct solution to original LP problem from bounded simplex output
            // Set non-basic variables to 0 or h based on result_init.e
            // Set basic variables to result_init.y0 or h-result_init.y0.

            // update DP_LPCA_problem.inB and DP_LPCA_problem.e by result_init.inB and result_init.e[0 to k=EffectorSize] that is e1(1:k+1) in matlab.  k+1 at all, so int (i=0;i<EffectorSize+1;++i) or (int i=0;i<DP_LPCA_problem.n;++i)
            for(int i=0;i<ControlSize;++i){
                DP_LPCA_problem.inB[i]=result_init.inB[i];
            }
            for(int i=0;i<DP_LPCA_problem.n;++i){
                DP_LPCA_problem.e[i]=result_init.e[i];
            }
            auto result = BoundedRevisedSimplex(DP_LPCA_problem);
            
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
                // std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
                // for (int i = 0; i < ControlSize; ++i) {
                //     std::cout << this->generalizedMoment[i] << std::endl;
                // }
            }
            if(result.errout)
            {
                err = 1;
                // std::cout << "Solver error"<< std::endl;
                // for (int i = 0; i < ControlSize; ++i) {
                //     std::cout << this->generalizedMoment[i] << std::endl;
                // }
            }
        }
        // Transform back to control variables
        for(int i=0;i<EffectorSize;++i){
            output[i]=xout[i]+this->aircraft.lowerLimits[i];
        }
        // Use upper_lam to prevent control surfaces from approaching position limits
        if(xout[EffectorSize]>1){
            for(int i=0;i<EffectorSize;++i){
                output[i]/=xout[EffectorSize];
            }
        }
        return;
    }
    
    
    void DPscaled_LPCA(float input[ControlSize], float output[EffectorSize], int& err, float & rho){
        // Direction Preserving Control Allocation Linear Program
        //     Reduced formulation (Solution Scaled from Boundary)

        // function [u,errout] = DPscaled_LPCA(yd,B,uMin,uMax,itlim);

        // Solves the control allocation problem while preserving the
        // objective direction for unattainable commands. The reduced
        // dimension of the linear program passed to the Bounded Revised
        // Simplex solver is formed by forcing the solution to be on the
        // boundary of the AMS and eliminating the highest magnitude
        // objective by solving the other constraints in terms of it.

        // For yd outside the AMS, the solution returned is that the
        // maximum in the direction of yd
        // B*u= lamda*yd
        // max lamda s.t. uMin <= u <= uMax

        // Reducing the degrees of freedom elminates the problems of redundant
        // solutions for attainable objectives. If the desired objective is on the
        // interior of the AMS the solution is scaled from the solution on the
        // boundary, yielding the same controls as the Direct Allocation solution.

        // (In the text this solution is discussed in section A.5.3)

        // (See Bodson, M., "Evaluation of Optimization Methods for
        //         Control Allocation",  AIAA 2001-4223).

        // Inputs:
        //         yd [n]    = Desired objective
        //         B [n,m]   = Control Effectiveness matrix
        //         uMin[m,1] = Lower bound for controls
        //         uMax[m,1] = Upper bound for controls
        //         itlim     = Number of allowed iterations limit
        //                         (Sum of iterations in both branches)

        // Outputs:
        //         u[m,1]     = Control Solution
        //         errout     = Error Status code
        //                         0 = found solution
        //                         <0 = Error in finding initial basic feasible solution
        //                         >0 = Error in finding final solution
        //                         -1,1 = Solver error (unbounded solution)
        //                         -2   = Initial feasible solution not found
        //                         -3,3 = Iteration limit exceeded

        // Calls:
        //         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)

        // Notes:
        // If yd is close to zero, u = 0;

        // Error code < 0 implies an error in the initialization and there is no guarantee on
        // the quality of the output solution other than the control limits.
        // Error code > 0 for errors in final solution.

        // Modification History
        // 2002      Roger Beck  Original ( DPcaLP2.m)
        // 8/2014    Roger Beck  Update
        // 4/2024    Meng ChaoHeng  Implement in cpp 
        // DPscaled_LPCA函数利用飞行器数据，将分配问题描述为DP_LP问题并用BoundedRevisedSimplex求解
        //=======================
        // std::cout << "DPscaled_LPCA"<< std::endl; 
        // Initialize error code to zero
        // int err = 0;
        // Figure out how big the problem is (use standard CA definitions for m & n)
        // but in here we use [m,k] = size(B) instead of [n,m] = size(B) in matlab. just for adapt to BoundedRevisedSimplex
        // we use [m,n] = size(A) in BoundedRevisedSimplex, that is, k + 1 = n.
        // Check to see if yd == 0
        // May want to adjust the tolerance to improve numerics of later steps
        bool flag=false;
        for(int i=0;i<ControlSize;++i){
            if(std::abs(input[i]) > DPscaled_LPCA_problem.tol)
            {
                flag = true; // Check this condition
                break;
            }
        }
        if(!flag){
            for(int i=0;i<EffectorSize;++i){
                output[i]=0;
            }
            // std::cout << "return 0"<< std::endl;
            return;
        }
        //=======================  
        // We inital this problem in constructor.
        // Construct an LP using scaling parameter to enforce direction preserving
        // To find Feasible solution construct problem with appended slack variables
        // ref. is A.6.4 Initialization of the Simplex Algorithm of <Aircraft control allocation>

        // now we update the problem by input data.
        // update A b c but not h, h = uMax-uMin is assumpt always the same.
        // float yd[3]={0.1,0.1,0.2}; // random value for inital.
        float my=input[0]; // 存储最大的绝对值
        int iy=0; // 存储最大绝对值的索引
        for (int i = 0; i < ControlSize; ++i) {
            float absValue = std::abs(input[i]); // 计算 yd 中第 i 个元素的绝对值
            if (absValue > my) { // 如果当前绝对值大于 my，则更新 my 和 iy
                my = absValue;
                iy = i;
            }
        }
        // std::cout << "DPscaled_LPCA_problem iy: "<< iy <<std::endl; 
        // std::cout << "DPscaled_LPCA_problem my: "<< my <<std::endl; 
        // copy firstly !!!
        float Bt[ControlSize][EffectorSize];
        float ydt[ControlSize];
        
        for(int i=0;i<ControlSize;++i){
            for (int j = 0; j <EffectorSize; ++j) {
                Bt[i][j] = this->aircraft.controlEffectMatrix[i][j];
            }
        }
        for(int i=0;i<ControlSize;++i){
            ydt[i]=input[i];
        }

        // move
        for(int j=0;j<EffectorSize;++j){
            Bt[0][j] = this->aircraft.controlEffectMatrix[iy][j];
            for (int i = iy; i >0; i--) {
                Bt[i][j] = this->aircraft.controlEffectMatrix[i-1][j];
            }
        }
        ydt[0] = input[iy];
        for (int j = iy; j >0; j--) {
            ydt[j] = input[j-1];
        }
        // swap 2 3 row of ydt
        float ydt2=ydt[1];
        ydt[1]=ydt[2];
        ydt[2]=ydt2;
        // swap 2 3 col of Bt
        for(int i=0;i<EffectorSize;++i){
            float Bt2=Bt[1][i];
            Bt[1][i]=Bt[2][i];
            Bt[2][i]=Bt2;
        }
        
        // std::cout << "aircraft.controlEffectMatrix:" << std::endl;
        // for (size_t i = 0; i < ControlSize; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << this->aircraft.controlEffectMatrix[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "Bt:" << std::endl;
        // for (size_t i = 0; i < ControlSize; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << Bt[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // M = [ydt(2:ControlSize) -ydt(1)*eye(ControlSize-1)];
        float M[ControlSize-1][ControlSize];
        // M[0][0]=ydt[1];
        // M[1][0]=ydt[2];
        // M[0][1]=-ydt[0];
        // M[1][1]=0;
        // M[0][2]=0;
        // M[1][2]=-ydt[0];
        // or 
        // 将 M 的所有元素初始化为0，并同时填充 M 的第一列和对角线元素
        for (int i = 0; i < ControlSize - 1; ++i) {
            for (int j = 0; j < ControlSize; ++j) {
                if (j == 0) {
                    // 填充 M 的第一列
                    M[i][0] = ydt[i + 1];
                } else if (j == i + 1) {
                    // 填充对角线元素
                    M[i][j] = -ydt[0];
                } else {
                    // 初始化其他元素为0
                    M[i][j] = 0.0f;
                }
            }
        }
        // std::cout << "M:" << std::endl;
        // for (size_t i = 0; i < ControlSize-1; ++i) {
        //     for (size_t j = 0; j < ControlSize; ++j) {
        //         std::cout << M[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        for (int i = 0; i < ControlSize-1; ++i) {
            for (int j = 0; j < EffectorSize; ++j) {
                DPscaled_LPCA_problem.A[i][j] = 0;
                for (int k = 0; k < ControlSize; ++k) {
                    DPscaled_LPCA_problem.A[i][j] += M[i][k] * Bt[k][j];
                }
            }
        }
        // std::cout << "DPscaled_LPCA_problem.A:" << std::endl;
        // for (size_t i = 0; i < ControlSize-1; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << DPscaled_LPCA_problem.A[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        for(int i=0; i<ControlSize-1; ++i)
        {
            float temp=0;
            for(int j=0; j<EffectorSize; ++j)
            {
                temp += -DPscaled_LPCA_problem.A[i][j]*this->aircraft.lowerLimits[j];
            }
            DPscaled_LPCA_problem.b[i] = temp;
        }
        for (int i = 0; i < EffectorSize; ++i) {
            DPscaled_LPCA_problem.c[i] = 0;
            for (int j = 0; j < ControlSize; ++j) {
                DPscaled_LPCA_problem.c[i] += -Bt[j][i] * ydt[j];
            }
        }
        for (int i = 0; i < EffectorSize; ++i) {
            DPscaled_LPCA_problem.h[i] = this->aircraft.upperLimits[i]-this->aircraft.lowerLimits[i];
        }
        // std::cout << "DPscaled_LPCA_problem b: [";
        // for (size_t i = 0; i < ControlSize-1; ++i) {
        //     std::cout << DPscaled_LPCA_problem.b[i];
        //     if (i < ControlSize-1 - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;


        // std::cout << "DPscaled_LPCA_problem c: [";
        // for (size_t i = 0; i <EffectorSize; ++i) {
        //     std::cout << DPscaled_LPCA_problem.c[i];
        //     if (i <EffectorSize - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        // std::cout << "DPscaled_LPCA_problem h: [";
        // for (size_t i = 0; i < EffectorSize; ++i) {
        //     std::cout << DPscaled_LPCA_problem.h[i];
        //     if (i < EffectorSize- 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        // update Ai and hi
        //==================================Pre_DPscaled_LPCA_problem================================
        // std::cout << "Pre_DPscaled_LPCA_problem"<< std::endl; 
        // Pre_DPscaled_LPCA_problem.tol=1e-7;
        // Pre_DPscaled_LPCA_problem.itlim = 10;
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            for(int j=0; j<DPscaled_LPCA_problem.n; ++j)
            {
                Pre_DPscaled_LPCA_problem.A[i][j] = DPscaled_LPCA_problem.A[i][j];
            }
            Pre_DPscaled_LPCA_problem.b[i] = DPscaled_LPCA_problem.b[i]; // the same as DPscaled_LPCA_problem

        }
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            Pre_DPscaled_LPCA_problem.A[i][i + DPscaled_LPCA_problem.n] = (DPscaled_LPCA_problem.b[i] > 0) ? 1 : -1; // sb = 2*(b > 0)-1; Ai = [A diag(sb)]; 
        }
        for(int i=0; i<DPscaled_LPCA_problem.n; ++i)
        {
            Pre_DPscaled_LPCA_problem.h[i] = DPscaled_LPCA_problem.h[i];
        }
        for(int i=0; i<DPscaled_LPCA_problem.m; ++i)
        {
            Pre_DPscaled_LPCA_problem.h[i+DPscaled_LPCA_problem.n] = 2*std::abs(DPscaled_LPCA_problem.b[i]);
        }
        // std::cout << "Pre_DPscaled_LPCA_problem A:" << std::endl;
        // for (size_t i = 0; i < Pre_DPscaled_LPCA_problem.m; ++i) {
        //     for (size_t j = 0; j < Pre_DPscaled_LPCA_problem.n; ++j) {
        //         std::cout << Pre_DPscaled_LPCA_problem.A[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "Pre_DPscaled_LPCA_problem b: [";
        // for (size_t i = 0; i < Pre_DPscaled_LPCA_problem.m; ++i) {
        //     std::cout << Pre_DPscaled_LPCA_problem.b[i];
        //     if (i < Pre_DPscaled_LPCA_problem.m - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        // std::cout << "Pre_DPscaled_LPCA_problem c: [";
        // for (size_t i = 0; i < Pre_DPscaled_LPCA_problem.n; ++i) {
        //     std::cout << Pre_DPscaled_LPCA_problem.c[i];
        //     if (i < Pre_DPscaled_LPCA_problem.n - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        // std::cout << "Pre_DPscaled_LPCA_problem h: [";
        // for (size_t i = 0; i < Pre_DPscaled_LPCA_problem.n; ++i) {
        //     std::cout << Pre_DPscaled_LPCA_problem.h[i];
        //     if (i < Pre_DPscaled_LPCA_problem.n- 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;
        // Use Bounded Revised Simplex to find initial basic feasible point of original program
        auto result_init = BoundedRevisedSimplex(Pre_DPscaled_LPCA_problem);

        // std::cout << "result_init"<< std::endl; 

        // Check that Feasible Solution was found
        if(result_init.iters>=Pre_DPscaled_LPCA_problem.itlim){
            err = -3;
            // std::cout << "Pre Too Many Iterations Finding Final Solution"<< std::endl; 
            // for (int i = 0; i < ControlSize; ++i) {
            //     std::cout << this->generalizedMoment[i] << std::endl;
            // }
        }
        for(int i=0;i<ControlSize-1;++i){
            if(result_init.inB[i]> EffectorSize-1) // DPscaled_LPCA_problem is origin problem, k=DPscaled_LPCA_problem.n-1 = EffectorSize 
            {
                // which mean inital basic index is out of the origin problem.
                err = -2;
                std::cout << "Pre No Initial Feasible Solution found"<< std::endl; 
                for (int k = 0; k < ControlSize; ++k) {
                    std::cout << this->generalizedMoment[k] << std::endl;
                }
                break;
            }
        }
        if(result_init.errout){
            err = -1;
            std::cout << "Pre Solver error"<< std::endl; 
            for (int i = 0; i < ControlSize; ++i) {
                std::cout << this->generalizedMoment[i] << std::endl;
            }
        }
        // solve Pre_DPscaled_LPCA_problem but proccess DPscaled_LPCA_problem
        float xout[EffectorSize];
        for(int i=0;i<DPscaled_LPCA_problem.n;++i){
            xout[i]=0;
        }
        if(err!=0) // Construct an incorrect solution to accompany error flags
        {
            // use result_init data
            // matlab: indv = inB1<=(k+1); xout(inB1(indv)) = y1(indv); % in matlab the index from 1 to k, but cpp is 0 to k-1
            for(int i=0;i<ControlSize-1;++i){
                xout[result_init.inB[i]]=result_init.y0[i];
                
            }
            // xout(~e1(1:k+1)) = -xout(~e1(1:k+1))+h(~e1(1:k+1));
            for(int i=0;i<EffectorSize;++i){
                if(!result_init.e[i]){
                    xout[i] = -xout[i] + DPscaled_LPCA_problem.h[i];
                }
            }
        }
        else //No Error continue to solve problem
        { 
            // Solve using initial problem from above
            // Construct solution to original LP problem from bounded simplex output
            // Set non-basic variables to 0 or h based on result_init.e
            // Set basic variables to result_init.y0 or h-result_init.y0.

            // update DPscaled_LPCA_problem.inB and DPscaled_LPCA_problem.e by result_init.inB and result_init.e[0 to k=EffectorSize] that is e1(1:k+1) in matlab.  k+1 at all, so int (i=0;i<EffectorSize+1;++i) or (int i=0;i<DPscaled_LPCA_problem.n;++i)
            for(int i=0;i<DPscaled_LPCA_problem.m;++i){
                DPscaled_LPCA_problem.inB[i]=result_init.inB[i];
            }
            for(int i=0;i<DPscaled_LPCA_problem.n;++i){
                DPscaled_LPCA_problem.e[i]=result_init.e[i];
            }
            // std::cout << "problem.inB: [";
            // for (size_t i = 0; i <DPscaled_LPCA_problem.m; ++i) {
            //     std::cout << DPscaled_LPCA_problem.inB[i];
            //     if (i <DPscaled_LPCA_problem.m- 1) {
            //         std::cout << ", ";
            //     }
            // }
            // std::cout << "]" << std::endl;
            // std::cout << "e: [";
            //     for (size_t i = 0; i < DPscaled_LPCA_problem.n; ++i) {
            //         std::cout << DPscaled_LPCA_problem.e[i];
            //         if (i < DPscaled_LPCA_problem.n- 1) {
            //             std::cout << ", ";
            //         }
            //     }
            // std::cout << "]" << std::endl;
            // std::cout << "before"<< std::endl; 
            auto result = BoundedRevisedSimplex(DPscaled_LPCA_problem);
            
            for(int i=0;i<ControlSize-1;++i){
                xout[result.inB[i]]=result.y0[i];
            }
            for(int i=0;i<EffectorSize;++i){
                if(!result.e[i]){
                    xout[i]=-xout[i]+DPscaled_LPCA_problem.h[i];
                }
            }

            if(result.iters>=DPscaled_LPCA_problem.itlim){
                err = 3;
                // std::cout << "Too Many Iterations Finding Final Solution"<< std::endl; 
                // for (int i = 0; i < ControlSize; ++i) {
                //     std::cout << this->generalizedMoment[i] << std::endl;
                // }
            }
            if(result.errout)
            {
                err = 1;
                // std::cout << "Solver error"<< std::endl;
                // for (int i = 0; i < ControlSize; ++i) {
                //     std::cout << this->generalizedMoment[i] << std::endl;
                // }
            }
        }
        // Transform back to control variables
        for(int i=0;i<EffectorSize;++i){
            output[i]=xout[i]+this->aircraft.lowerLimits[i];
            // std::cout << "output[i]"<<output[i]<< std::endl;
        }
        // Use upper_lam to prevent control surfaces from approaching position limits
        rho = calculateRho(ydt, output, Bt, DPscaled_LPCA_problem.tol);
        // std::cout << "Bt:" << std::endl;
        // for (size_t i = 0; i < ControlSize; ++i) {
        //     for (size_t j = 0; j < EffectorSize; ++j) {
        //         std::cout << Bt[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << "ydt: [";
        // for (size_t i = 0; i < ControlSize; ++i) {
        //     std::cout << ydt[i];
        //     if (i < ControlSize - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << "]" << std::endl;

        // std::cout << "rho"<<rho<< std::endl;
        if(rho>1){
            for(int i=0;i<EffectorSize;++i){
                output[i]/=rho;
            }
        }
        return;
    }
    
    // 其他成员函数和成员变量定义
    float generalizedMoment[ControlSize]; // 构造函数设置
    // 线性规划相关
    LinearProgrammingProblem<ControlSize, EffectorSize+1> DP_LPCA_problem;// 提前设置 inital by  aircraft data 
    LinearProgrammingProblem<ControlSize, (EffectorSize+1) + ControlSize> Pre_DP_LPCA_problem;// 提前设置 inital by aircraft data

    LinearProgrammingProblem<ControlSize-1, EffectorSize> DPscaled_LPCA_problem;// 提前设置 inital by  aircraft data 
    LinearProgrammingProblem<ControlSize-1, EffectorSize + (ControlSize-1)> Pre_DPscaled_LPCA_problem;// 提前设置 inital by aircraft data
    float upper_lam = 1e4;
    
};
//==================
// DP_LPCA.allocateControl_bases_solver execution time: 6.67e-07s
// DP_LPCA.allocateControl execution time: 7.5e-07s
//==================

// and user can define more...
// template <int ControlSize, int EffectSize>
// class XX_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectSize> {}