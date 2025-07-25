cmake_minimum_required(VERSION 3.0)
project(MyProject)

# 设置编译器标准
set(CMAKE_CXX_STANDARD 14)# for matrix
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置编译器选项（包括-O3优化）
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# 原生 ARM64 编译​（移除 -march=native 等标志）
# 添加编译器标志以启用SIMD指令集
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
# ✅ 安全：仅 Apple 平台生效
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "arm64")
endif()

# 添加头文件搜索路径
include_directories(${CMAKE_SOURCE_DIR}/alglib)
include_directories(${CMAKE_SOURCE_DIR}/src/ControlAllocation) # 如果您有其他的头文件目录
include_directories(${CMAKE_SOURCE_DIR}/src/LPwrap) # 如果您有其他的头文件目录
include_directories(${CMAKE_SOURCE_DIR}/src/QP) # 如果您有其他的头文件目录
# 查找Eigen库并包含它
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
# or
# include_directories(/usr/local/Cellar/eigen/3.4.0_1/include/eigen3)

# 添加子目录alglib
add_subdirectory(alglib)

# 添加 matrix 库的子目录
# add_subdirectory(${CMAKE_SOURCE_DIR}/matrix) # maybe don't need

# 将 matrix 库的头文件路径添加到包含目录
include_directories(${CMAKE_SOURCE_DIR}/matrix)


# 添加可执行文件
add_executable(main test/main.cpp src/LPwrap/allocator_dir_LPwrap_4.c src/LPwrap/rt_nonfinite.c src/QP/wls_alloc_gen.c src/QP/wls_alloc_gen_emxutil.c)
add_executable(alglib_based_minlp_basic test/alglib_based_minlp_basic.cpp)
add_executable(matlab_code_gen_LPwrap_test test/matlab_code_gen_LPwrap_test.cpp src/LPwrap/allocator_dir_LPwrap_4.c src/LPwrap/rt_nonfinite.c)
add_executable(Eingen_based_simplex test/Eingen_based_simplex.cpp src/LPwrap/allocator_dir_LPwrap_4.c src/LPwrap/rt_nonfinite.c)
# 添加链接的库
target_link_libraries(main alglib)
target_link_libraries(alglib_based_minlp_basic alglib)
target_link_libraries(matlab_code_gen_LPwrap_test alglib)
target_link_libraries(Eingen_based_simplex  alglib)
