# About
 This project refers to the `control_allocation_lib` of [control allocation library](https://github.com/mengchaoheng/control_allocation.git) , especially the code of [ simulation of the book "aircraft control allocation"](https://github.com/mengchaoheng/aircraft-control-allocation), and using C++ to implement the LP-based control allocation algorithm.
 ## Dependencies
 This project depends on [PX4-matrix](https://github.com/mengchaoheng/PX4-Matrix.git). PX4-matrix is a submodule of this project now. If you want to apply these algorithms to your project, please clone the PX4-matrix library to your project first, and copy `alloc_based_LP_cpp/src/ControlAllocation/ControlAllocation.h` to the project , just include ControlAllocation.h when using it.
### download alglib (option)
If you want to test the algorithm based alglib, download and Unzip alglib to the root folder of this project. rename as `alglib`. 

## Build and run
```sh
cd build
cmake ..
make 
```

cd to build path, and run
```sh
./main ## or other target: matlab_code_gen_LPwrap_test, Eingen_based_simplex, alglib_based_minlp_basic.
```
then the main target will write data to `output.csv`, and you can run the `test.m` to compare the output to that generate by matlab code.
## Usage
1. Initialize an aircraft object and allocator object.
```C++
Aircraft<3, 4> df_4(_B, -0.3491, 0.3491); // Create an aircraft object with 4 steering vectors and 3 generalized moments
DP_LP_ControlAllocator<3, 4> Allocator(df_4); // Create a control distributor object for an aircraft with 4 steering vectors and 3 generalized moments (translated into a linear programming problem with dimensions related to parameters <3, 4>.)
```
2. After the input is updated, run `Allocator.DP_LPCA(input, u, err)` for example.

## Test report
running time test on MacOS

Eigen based:
```
DPscaled_LPCA execution time: 2.0875e-05s
DP_LPCA execution time: 3.5875e-05s
```

the Matrix lib of px4 based:
```
Allocator.allocateControl Average execution time: 7.1138e-07s
Allocator.DPscaled_LPCA Average execution time: 5.03832e-07s
Allocator.DP_LPCA Average execution time: 1.50543e-06s
allocator_dir_LPwrap_4 Average execution time: 1.41764e-06s
```