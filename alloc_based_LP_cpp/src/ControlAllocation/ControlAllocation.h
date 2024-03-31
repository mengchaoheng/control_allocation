
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

// template<int M, int N>
// class LinearProgrammingSolver {
// public:
//     LinearProgrammingSolver() = default;
//     // 构造函数，接受一个 LinearProgrammingProblem 对象作为参数
//     LinearProgrammingSolver(const LinearProgrammingProblem<M, N>& inputProblem) : problem(inputProblem) {}

//     // 成员函数调用 BoundedRevisedSimplex 算法，并返回结果
//     void solve() {
//         result = BoundedRevisedSimplex(problem);
//     }
//     LinearProgrammingProblem<M, N> problem;
//     LinearProgrammingResult<M, N> result;

//     // 获取结果
//     const LinearProgrammingResult<M, N>& getResult() const {
//         return result;
//     }
// };


// 定义函数模板
template<int M, int N>
LinearProgrammingResult<M, N> BoundedRevisedSimplex(LinearProgrammingProblem<M, N>& problem) {
    LinearProgrammingResult<M, N> result;
    // 实现线性规划算法
    // 使用 problem.inB, problem.inD, problem.itlim, problem.A, problem.b, problem.c, problem.h, problem.e
    std::cout << "Solver error"<< std::endl;
    float tol=1e-7;
    const int n_m=N-M;
    int* nind = problem.generateSequence(0, n_m-1);
    
    int* ind_all = problem.generateSequence(0, N-1);
    problem.setdiff(ind_all, N, problem.inB, M, problem.inD);
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




// 飞行器基类模板
template <int EffectSize, int ControlSize>
class AircraftBase {
public:
    float controlVector[ControlSize]; // 操纵向量
    float controlEffectMatrix[EffectSize][ControlSize]; // 控制效应矩阵 (generalizedMomentSize X controlVectorSize)
    float generalizedMoment[EffectSize]; // 广义力矩
    float upperLimits[ControlSize]; // 操纵向量上限变量
    float lowerLimits[ControlSize]; // 操纵向量下限变量

    // 构造函数
    // 拷贝构造函数
    AircraftBase(const AircraftBase& other) {
        // 复制 controlVector
        for (int i = 0; i < ControlSize; ++i) {
            controlVector[i] = other.controlVector[i];
            upperLimits[i] = other.upperLimits[i];
            lowerLimits[i] = other.lowerLimits[i];
        }

        // 复制 controlEffectMatrix
        for (int i = 0; i < EffectSize; ++i) {
            for (int j = 0; j < ControlSize; ++j) {
                controlEffectMatrix[i][j] = other.controlEffectMatrix[i][j];
            }
        }

        // 复制 generalizedMoment
        for (int i = 0; i < EffectSize; ++i) {
            generalizedMoment[i] = other.generalizedMoment[i];
        }
    }


    AircraftBase() {
        // 初始化 controlVector, controlEffectMatrix, generalizedMoment, upperLimits, lowerLimits 等数组
        // 可以使用默认初始化或者自定义初始化方式
        // 例如：
        for (int i = 0; i < ControlSize; ++i) {
            controlVector[i] = 0.0f;
            upperLimits[i] = 0.0f;
            lowerLimits[i] = 0.0f;
            for (int j = 0; j < EffectSize; ++j) {
                controlEffectMatrix[j][i] = 0.0f;
            }
        }
        for (int i = 0; i < EffectSize; ++i) {
            generalizedMoment[i] = 0.0f;
        }
    }

    // 析构函数
    ~AircraftBase() {
        // 不需要手动释放内存，因为数组是在栈上分配的，会在对象生命周期结束时自动释放
    }
};

// 飞行器类模板, 不同飞机定义新的类继承基类
template <int EffectSize, int ControlSize>
class Aircraft : public AircraftBase<ControlSize, EffectSize> {
private:
    // 添加特定飞行器类型的模型参数
public:
    // 构造函数
    // 拷贝构造函数
    // 拷贝构造函数
    Aircraft(const Aircraft& other) : AircraftBase<ControlSize, EffectSize>(other) {
        // 将其他对象的成员变量值复制到新对象中
        l1 = other.l1;
        l2 = other.l2;
        k_v = other.k_v;
        upper = other.upper;
        lower = other.lower;
    }
    Aircraft() : AircraftBase<ControlSize, EffectSize>() {
        // 可选的初始化代码
        l1=0;
        l2=0;
        k_v=0;
    }
    // 构造函数，接受对应于飞行器类模板参数的初始化参数
    Aircraft(float controlEffectMatrixInit[EffectSize][ControlSize], 
             float upperLimitsInit[ControlSize], 
             float lowerLimitsInit[ControlSize]) {
        // 使用传入的初始化参数对飞行器的数组成员进行初始化
        for (int i = 0; i < ControlSize; ++i) {
            this->controlVector[i] = 0;
            this->upperLimits[i] = upperLimitsInit[i];
            this->lowerLimits[i] = lowerLimitsInit[i];
            for (int j = 0; j < EffectSize; ++j) {
                this->controlEffectMatrix[j][i] = controlEffectMatrixInit[j][i];
            }
        }
        for (int i = 0; i < EffectSize; ++i) {
            this->generalizedMoment[i] = 0;
        }
    }
    // 构造函数，接受对应于飞行器类模板参数的初始化参数
    Aircraft(float upperLimitsInit[ControlSize], 
             float lowerLimitsInit[ControlSize]) {
        // 使用传入的初始化参数和飞行器
        for (int i = 0; i < ControlSize; ++i) {
            this->controlVector[i] = 0;
            this->upperLimits[i] = upperLimitsInit[i];
            this->lowerLimits[i] = lowerLimitsInit[i];
            for (int j = 0; j < EffectSize; ++j) {
                this->controlEffectMatrix[j][i] = 0;
            }
        }
        for (int i = 0; i < EffectSize; ++i) {
            this->generalizedMoment[i] = 0;
        }
        // and define other value manual to set B (controlEffectMatrix).
    }
    Aircraft(float controlEffectMatrixInit[EffectSize][ControlSize], float upper, float lower) :  upper(upper), lower(lower){
        // 使用传入的初始化参数和飞行器
        for (int i = 0; i < ControlSize; ++i) {
            this->controlVector[i] = 0;
            this->upperLimits[i] = upper;
            this->lowerLimits[i] = lower;
            for (int j = 0; j < EffectSize; ++j) {
                this->controlEffectMatrix[j][i] = controlEffectMatrixInit[j][i];
            }
        }
        for (int i = 0; i < EffectSize; ++i) {
            this->generalizedMoment[i] = 0;
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
template <int ControlSize, int EffectSize>
class ControlAllocatorBase {
public:
    // 默认构造函数
    ControlAllocatorBase() {
        // 在此初始化成员变量，或者留空
    }
    // 参数列表构造函数
    ControlAllocatorBase(const Aircraft<ControlSize, EffectSize>& aircraft, float generalizedMoment[ControlSize]) {
        // 在此初始化成员变量，可以使用传入的参数值
        this->aircraft = aircraft;
        for (int i = 0; i < ControlSize; ++i) {
            this->generalizedMoment[i] = generalizedMoment[i];
        }
    }

    virtual float*  allocateControl(float generalizedMoment[ControlSize]) = 0;

    // 其他数学函数和成员变量定义
    Aircraft<ControlSize, EffectSize> aircraft;// 构造函数设置
    float generalizedMoment[ControlSize];// 构造函数设置
};
// 控制分配类模板
template <int ControlSize, int EffectSize>
class DP_LP_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectSize> {
private:
    // 添加算法设置参数
public:
    // 构造函数, 利用aircraft 预设置 LinearProgrammingSolver 和 LinearProgrammingProblem
    // 构造函数
    DP_LP_ControlAllocator(const Aircraft<ControlSize, EffectSize>& aircraft, float generalizedMoment[ControlSize])
        : ControlAllocatorBase<ControlSize, EffectSize>(aircraft, generalizedMoment){
        // 在此处初始化其他成员变量 DP_LPCA_problem 和 Pre_DP_LPCA_problem 
        // 线性规划数据
        DP_LPCA_problem.itlim = 100;
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
            DP_LPCA_problem.A[i][DP_LPCA_problem.n-1] =-generalizedMoment[i]; // will be change in the loop
            this->aircraft.generalizedMoment[i]=generalizedMoment[i];
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
            // 
            upper_lam=1e4;
        }

    // 设置算法参数函数

    // 析构函数
    ~DP_LP_ControlAllocator() {
        // 如果有需要释放的资源，可以在这里添加代码
    }

    float* allocateControl(float generalizedMoment[ControlSize]) override {
        // 重写控制分配器函数
        // 实现控制分配算法
        float* control = new float[ControlSize];
        // 算法实现
        // DP_LPCA（generalizedMoment, aircraft） 
        // DP_LPCA函数利用飞行器数据，将分配问题描述为DP_LP问题并用BoundedRevisedSimplex求解
        // 使用模版函数result = BoundedRevisedSimplex(problem);
        this->DP_LPCA_problem.A[0][DP_LPCA_problem.n-1]=-generalizedMoment[0];
        this->DP_LPCA_problem.A[1][DP_LPCA_problem.n-1]=-generalizedMoment[1];
        this->DP_LPCA_problem.A[2][DP_LPCA_problem.n-1]=-generalizedMoment[2];
        auto result = BoundedRevisedSimplex(this->DP_LPCA_problem);
        // for(int i==0;i<ControlSize;++i)
        // {
        //     control[i]=result
        // }
        // 使用结果
        // result.y0, result.inB, result.e, result.errout
        int err = 0;
        float xout[DP_LPCA_problem.n];
        for(int i=0;i<DP_LPCA_problem.n;++i){
            xout[i]=0;
        }
        for(int i=0;i<DP_LPCA_problem.m;++i){
            xout[result.inB[i]]=result.y0[i];
        }
        std::cout << "result.y0: [";
    for (size_t i = 0; i < DP_LPCA_problem.n-1; ++i) {
        std::cout << result.y0[i];
        if (i < DP_LPCA_problem.n-2) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

        for(int i=0;i<DP_LPCA_problem.n;++i){
            if(!result.e[i]){
                xout[i]=-xout[i]+DP_LPCA_problem.h[i];
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
        for(int i=0;i<ControlSize;++i){
            control[i]=xout[i]+this->aircraft.lowerLimits[i];
        }
        if(xout[ControlSize]>1){
            for(int i=0;i<ControlSize;++i){
                control[i]/=xout[ControlSize];
            }
        }

        return control;
    }

    // 其他成员函数和成员变量定义
    LinearProgrammingProblem<ControlSize, EffectSize+1> DP_LPCA_problem;// 提前设置 using aircraft
    LinearProgrammingProblem<ControlSize, EffectSize+4> Pre_DP_LPCA_problem;// 提前设置 using aircraft

    float upper_lam;


    // LinearProgrammingSolver<ControlSize, EffectSize+1> LPSolver;// 提前设置 using DP_LPCA_problem
    // LinearProgrammingSolver<ControlSize, EffectSize+1> Pre_LPSolver;// 提前设置 using Pre_DP_LPCA_problem

    // void DP_LPCA(float generalizedMoment[ControlSize], const AircraftBase<ControlSize, EffectSize>& aircraft);//using LPSolver and Pre_LPSolver

};

// template <int ControlSize, int EffectSize>
// class XX_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectSize> {}