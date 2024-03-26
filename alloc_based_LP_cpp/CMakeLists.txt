cmake_minimum_required(VERSION 3.0)
project(MyProject)

# 设置编译器标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加头文件搜索路径
include_directories(${CMAKE_SOURCE_DIR}/alglib)
include_directories(${CMAKE_SOURCE_DIR}/src/LPWrap) # 如果您有其他的头文件目录

# 查找Eigen库并包含它
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# 添加子目录alglib
add_subdirectory(alglib)

# 添加可执行文件
add_executable(MyExecutable src/main.cpp)
add_executable(minlp_basic src/minlp_basic.cpp)
add_executable(LPwrap_test src/LPWrap/LPwrap_test.cpp src/LPWrap/allocator_dir_LPwrap_4.c src/LPWrap/rt_nonfinite.c)

# 添加链接的库
target_link_libraries(MyExecutable Eigen3::Eigen alglib)
target_link_libraries(minlp_basic alglib)
target_link_libraries(LPwrap_test alglib)
