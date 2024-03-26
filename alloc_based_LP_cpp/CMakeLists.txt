cmake_minimum_required(VERSION 3.0)
project(EigenExample)

# 设置编译器标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 寻找Eigen库的路径
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# 添加可执行文件
add_executable(main main.cpp)
