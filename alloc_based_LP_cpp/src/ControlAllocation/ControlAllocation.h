
// #include <matrix/math.hpp>
// using namespace matrix;
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
    ControlAllocatorBase(Aircraft<ControlSize, EffectSize>& aircraft, float generalizedMoment[ControlSize]) {
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
        // 你可以在这里为它们设置其他值
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
        return control;
    }

    // 其他成员函数和成员变量定义
    LinearProgrammingProblem<ControlSize, EffectSize+1> DP_LPCA_problem;// 提前设置 using aircraft
    LinearProgrammingProblem<ControlSize, EffectSize+4> Pre_DP_LPCA_problem;// 提前设置 using aircraft

    // LinearProgrammingSolver<ControlSize, EffectSize+1> LPSolver;// 提前设置 using DP_LPCA_problem
    // LinearProgrammingSolver<ControlSize, EffectSize+1> Pre_LPSolver;// 提前设置 using Pre_DP_LPCA_problem

    void DP_LPCA(float generalizedMoment[ControlSize], const AircraftBase<ControlSize, EffectSize>& aircraft);//using LPSolver and Pre_LPSolver

};
