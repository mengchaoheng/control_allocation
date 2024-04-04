## download alglib
Unzip to the root folder of this project. rename as `alglib`.

## build and run
```sh
cd build
cmake ..
make && ./main ## or other target: matlab_code_gen_LPwrap_test, Eingen_based_simplex, alglib_based_minlp_basic.
```

cd to build path, and run
```sh
./main ## or other target: matlab_code_gen_LPwrap_test, Eingen_based_simplex, alglib_based_minlp_basic.
```

## test report
running time test on MacOS

Eigen based:
```
DPscaled_LPCA execution time: 2.0875e-05s
DP_LPCA execution time: 3.5875e-05s
```

the Matrix lib of px4 based:
```
allocator_dir_LPwrap_4 execution time: 6.67e-07s
Allocator.allocateControl execution time: 1.66e-07s
Allocator.allocateControl_bases_solver execution time: 6.66e-07s
Allocator.DPscaled_LPCA execution time: 4.17e-07s
Allocator.DP_LPCA execution time: 1.583e-06s
```