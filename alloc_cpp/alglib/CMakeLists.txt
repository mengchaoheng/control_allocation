# 添加源文件
add_library(alglib
    src/alglibinternal.cpp   
    src/dataanalysis.cpp     
    src/integration.cpp      
    src/kernels_fma.cpp      
    src/optimization.cpp     
    src/statistics.cpp
    src/alglibmisc.cpp       
    src/diffequations.cpp    
    src/interpolation.cpp    
    src/kernels_sse2.cpp     
    src/solvers.cpp
    src/ap.cpp               
    src/fasttransforms.cpp   
    src/kernels_avx2.cpp     
    src/linalg.cpp           
    src/specialfunctions.cpp
)


# 将include文件夹添加到头文件搜索路径
target_include_directories(alglib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
