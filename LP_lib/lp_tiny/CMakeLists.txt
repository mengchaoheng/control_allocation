cmake_minimum_required(VERSION 3.0)
project(test_input_alloc)
 

# 尝试自动发现LAPACK和BLAS包
find_package(LAPACK)
find_package(BLAS)
 
# 如果找到了LAPACK和BLAS，将它们链接到你的目标
if(LAPACK_FOUND AND BLAS_FOUND)
    add_executable(test_input_alloc test_input_alloc.cpp lp_tiny.c)
    target_link_libraries(test_input_alloc ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
else()
    message(FATAL_ERROR "LAPACK and/or BLAS not found.")
endif()