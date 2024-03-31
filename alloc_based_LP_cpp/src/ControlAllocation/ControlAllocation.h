
// #include <matrix/math.hpp>
// using namespace matrix;

// 飞行器基类模板
template <int ControlSize, int EffectSize>
class AircraftBase {
public:
    float controlVector[ControlSize]; // 操纵向量
    float controlEffectMatrix[EffectSize][ControlSize]; // 控制效应矩阵 (generalizedMomentSize X controlVectorSize)
    float generalizedMoment[EffectSize]; // 广义力矩
    float upperLimits[ControlSize]; // 操纵向量上限变量
    float lowerLimits[ControlSize]; // 操纵向量下限变量

    // 构造函数
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
template <int ControlSize, int EffectSize>
class Aircraft : public AircraftBase<ControlSize, EffectSize> {
private:
    // 添加特定飞行器类型的模型参数
public:
    // 构造函数

    // 设置模型参数函数

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
    virtual float*  allocateControl(const AircraftBase<ControlSize, EffectSize>& aircraft) = 0;

    // 其他数学函数和成员变量定义
};
// 控制分配类模板
template <int ControlSize, int EffectSize>
class ControlAllocator : public ControlAllocatorBase<ControlSize, EffectSize> {
private:
    // 添加算法设置参数
    // 添加飞行器类成员变量
public:
    // 构造函数

    // 设置算法参数函数

    // 析构函数
    ~ControlAllocator() {
        // 如果有需要释放的资源，可以在这里添加代码
    }

    float* allocateControl(const AircraftBase<ControlSize, EffectSize>& aircraft) override {
        // 重写控制分配器函数
        // 实现控制分配算法
        float* control = new float[ControlSize];
        // 算法实现
        return control;
    }

    // 其他成员函数和成员变量定义
};
