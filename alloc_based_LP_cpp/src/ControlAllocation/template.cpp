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
    int itlim;
    bool errout;
    // 其他结果成员
    // 默认构造函数，将所有成员变量初始化为0
    LinearProgrammingResult() : itlim(0), errout(false) {
        // 将数组成员变量初始化为0
        for (int i = 0; i < M; ++i) {
            y0[i] = 0.0f;
            inB[i] = 0;
        }
        for (int i = 0; i < N; ++i) {
            e[i] = false;
        }
    }
};




template<int M, int N>
class LinearProgrammingSolverForAC {
public:
    // 默认构造函数，初始化所有成员变量
    LinearProgrammingSolverForAC() 
        : problem(),                // 使用默认构造函数初始化 problem
          result(),                 // 使用默认构造函数初始化 result
          A_inB(),                  // 使用默认构造函数初始化 A_inB
          A_inD(),                  // 使用默认构造函数初始化 A_inD
          c_inB(),                  // 使用默认构造函数初始化 c_inB
          c_inD(),                  // 使用默认构造函数初始化 c_inD
          lamt(),                   // 使用默认构造函数初始化 lamt
          rdt(),                    // 使用默认构造函数初始化 rdt
          A_qel(),                  // 使用默认构造函数初始化 A_qel
          yq(),                     // 使用默认构造函数初始化 yq
          rat(),                    // 使用默认构造函数初始化 rat
          tol(1e-7),                // 将 tol 初始化为 1e-7
          n_m(N - M),               // 计算并初始化 n_m
          nind(generateSequence(0, N - M - 1)), // 使用 generateSequence 初始化 nind
          ind_all(generateSequence(0, N - 1))    // 使用 generateSequence 初始化 ind_all
    {
        // 构造函数体为空，所有成员变量已在初始化列表中初始化
    }
    // 构造函数，接受一个 LinearProgrammingProblem 对象作为参数
    LinearProgrammingSolverForAC(const LinearProgrammingProblem<M, N>& inputProblem) : problem(inputProblem) {
    }
    // 成员函数调用 Simplex （改自 BoundedRevisedSimplex ），并返回结果
    void solve() {
        result = Simplex(problem);
    }
    LinearProgrammingProblem<M, N> problem;
    LinearProgrammingResult<M, N> result;
    // simplex
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
    int* nind = generateSequence(0, N-M-1);
    int* ind_all = generateSequence(0, N-1);
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
    ControlAllocatorBase(const Aircraft<ControlSize, EffectorSize>& aircraft, float generalizedMoment[ControlSize]) 
        : aircraft(aircraft) { // 使用传入的aircraft对象初始化aircraft成员
        // 在此初始化成员变量，可以使用传入的参数值
        for (int i = 0; i < ControlSize; ++i) {
            this->generalizedMoment[i] = generalizedMoment[i];
        }
    }

    virtual float*  allocateControl(float generalizedMoment[ControlSize]) = 0;

    // 其他数学函数和成员变量定义
    Aircraft<ControlSize, EffectorSize> aircraft; // 构造函数设置
    float generalizedMoment[ControlSize]; // 构造函数设置
};

template <int ControlSize, int EffectorSize>
class DP_LP_ControlAllocator : public ControlAllocatorBase<ControlSize, EffectorSize> {
private:
    // 添加算法设置参数
public:
    // 构造函数, 利用aircraft 预设置LinearProgrammingProblem, 再初始化 LinearProgrammingSolverForAC
    // 构造函数
    DP_LP_ControlAllocator(const Aircraft<ControlSize, EffectorSize>& aircraft, float generalizedMoment[ControlSize])
        : ControlAllocatorBase<ControlSize, EffectorSize>(aircraft, generalizedMoment){
        // 在此处用aircraft, generalizedMoment初始化 成员变量 DP_LPCA_problem 和 Pre_DP_LPCA_problem 
        // 线性规划数据
        
        
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
        return EffectorCommand;
    }
    float* allocateControl_bases_solver(float input[ControlSize]){

    }

    // 其他成员函数和成员变量定义
    // 线性规划相关
    LinearProgrammingProblem<ControlSize, EffectorSize+1> DP_LPCA_problem;// 提前设置 using aircraft
    LinearProgrammingProblem<ControlSize, EffectorSize+4> Pre_DP_LPCA_problem;// 提前设置 using aircraft

    LinearProgrammingSolverForAC<ControlSize, EffectorSize+1> LPsolverForAC;

    float upper_lam = 1e4;
};
DP_LP_ControlAllocator的构造函数设置了DP_LPCA_problem的部分成员值，它的其他成员如何初始化的？