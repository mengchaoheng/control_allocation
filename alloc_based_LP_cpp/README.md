## download alglib
Unzip to the root folder of this project. rename as `alglib`.

## build and run
```sh
cd build
cmake ..
make && ./new_simplex ## or other target
```

cd to build path, and run
```sh
./new_simplex ## or other target
```

## test report
running time test on MacOS
```
DPscaled_LPCA execution time: 2.0875e-05s
DP_LPCA execution time: 3.5875e-05s
allocator_dir_LPwrap_4 execution time: 6.67e-07s
```