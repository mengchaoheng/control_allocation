# About
 This project refers to the `control_allocation_lib` of [control allocation library](https://github.com/mengchaoheng/control_allocation.git) , especially the code of [ simulation of the book "aircraft control allocation"](https://github.com/mengchaoheng/aircraft-control-allocation), and using C++ to implement the LP-based control allocation algorithm.
 ## Dependencies
 This project depends on [PX4-matrix](https://github.com/mengchaoheng/PX4-Matrix.git). PX4-matrix is a submodule of this project now. If you want to apply these algorithms to your project, please clone the PX4-matrix library to your project first, and copy `alloc_cpp/src/ControlAllocation/ControlAllocation.h` to the project , just include ControlAllocation.h when using it.
### download alglib (option)
If you want to test the algorithm based alglib, download and Unzip alglib to the root folder of this project. rename as `alglib`. 

## Build and run
To get input data:

1. Run `plot_fly_log_states.m` to get the `flight.mat` file.

2. Run `Generate_input_data.m` t get input data `input.csv`.
 and then 
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
2. After the input is updated, execute `Allocator.DP_LPCA(input, u, err)`, or `Allocator.allocateControl`, `Allocator.DPscaled_LPCA` for example.

## Test report

### Running time test on MacOS

1. Eigen based:
```
DPscaled_LPCA execution time: 2.0875e-05s
DP_LPCA execution time: 3.5875e-05s
```

2. The Matrix lib of px4 based:
```sh
# DP_LPCA and DP_LPCA_prio DPscaled_LPCA run restoring
Allocator.allocateControl Average execution time: 4.95833e-07s
Allocator.DPscaled_LPCA Average execution time: 7.42731e-07s
Allocator.DP_LPCA Average execution time: 1.34781e-06s
allocator_dir_LPwrap_4 Average execution time: 8.90307e-07s
Allocator.DP_LPCA_prio Average execution time: 1.34477e-06s

# DP_LPCA and DP_LPCA_prio DPscaled_LPCA  remove restoring
Allocator.DPscaled_LPCA Average execution time: 6.29866e-07s
Allocator.DP_LPCA Average execution time: 1.12595e-06s
Allocator.DP_LPCA_prio Average execution time: 1.11555e-06s

# add a higher unattainable for prio under restoring
Allocator.DP_LPCA_prio Average execution time: 2.01469e-06s
```
Update: add an update fuction, test result:
```sh
Allocator.allocateControl Average execution time: 4.22216e-07s
Allocator.DPscaled_LPCA Average execution time: 6.1186e-07s
Allocator.DP_LPCA Average execution time: 1.22761e-06s
allocator_dir_LPwrap_4 Average execution time: 8.1512e-07s
Allocator.DP_LPCA_prio Average execution time: 1.23046e-06s
wls_alloc_gen Average execution time: 5.63349e-07s
```

### Running on pixhawk

allocate the input[3], the running time is (us):
```c++
// =====================run test for allocation running time===========================================
	float m_higher[3]={0.0,  0.0,  0.0f}; //
	float m_lower[3]={20.0f,  10.0f,   0.0f};
	float input[3]={m_lower[0]+m_higher[0],  m_lower[1]+m_higher[1],   m_lower[2]+m_higher[2]};
	//==========================allocateControl===========================
	float u1[4]; int err1=0;
	timestamp_ca_start = hrt_absolute_time();
	Allocator.allocateControl(input, u1, err1);
	timestamp_ca_end = hrt_absolute_time();
	PX4_INFO("allocateControl: u1: %f, u2: %f, u3: %f, u4: %f. \n",(double) u1[0],(double) u1[1],(double) u1[2],(double) u1[3]);
	PX4_INFO("allocateControl time: %lld \n", (timestamp_ca_end - timestamp_ca_start) ); //nuttx
	//=========================DPscaled_LPCA============================INFO  [mixer_module] dir_alloc_sim time: 16
	float u2[4];int err2=0;float rho2=0;float u2_tmp[4];
	timestamp_ca_start = hrt_absolute_time();
	Allocator.DPscaled_LPCA(input, u2_tmp, err2, rho2);
	Allocator.restoring(u2_tmp,u2);
	timestamp_ca_end = hrt_absolute_time();
	PX4_INFO("DPscaled_LPCA: u1: %f, u2: %f, u3: %f, u4: %f. \n",(double) u2[0],(double) u2[1],(double) u2[2],(double) u2[3]);
	PX4_INFO("DPscaled_LPCA time: %lld \n", (timestamp_ca_end - timestamp_ca_start) ); //nuttx
	//========================DP_LPCA=============================
	float u3[4];int err3=0;float rho3=0;float u3_tmp[4];
	timestamp_ca_start = hrt_absolute_time();
	Allocator.DP_LPCA(input, u3_tmp, err3, rho3);
	Allocator.restoring(u3_tmp,u3);
	timestamp_ca_end = hrt_absolute_time();
	PX4_INFO("DP_LPCA: u1: %f, u2: %f, u3: %f, u4: %f. \n",(double) u3[0],(double) u3[1],(double) u3[2],(double) u3[3]);
	PX4_INFO("DP_LPCA time: %lld \n", (timestamp_ca_end - timestamp_ca_start) ); //nuttx
	//========================DP_LPCA_prio=============================
	float u4[4];int err4=0;float rho4=0;float u4_tmp[4];
	// float m_lower[3]={30.0f,  0.0f,   -0.0f};
        float u5[4];int err5=0;float rho5=0;float u5_tmp[4]; float m_tmp[3]={0.0,  0.0,  0.0f};
	timestamp_ca_start = hrt_absolute_time();
        Allocator.DP_LPCA_copy(m_higher,m_lower, u4_tmp, err4, rho4);
        if (err4<0){
            Allocator.DP_LPCA_copy(m_tmp,m_higher, u5_tmp, err5, rho5);
            Allocator.restoring(u5_tmp,u4);
        }else{
            Allocator.restoring(u4_tmp,u4);
        }
	timestamp_ca_end = hrt_absolute_time();
	PX4_INFO("DP_LPCA_prio: u1: %f, u2: %f, u3: %f, u4: %f. \n",(double) u4[0],(double) u4[1],(double) u4[2],(double) u4[3]);
	PX4_INFO("DP_LPCA_prio time: %lld \n", (timestamp_ca_end - timestamp_ca_start) ); //nuttx
```
without restoring
```sh 
INFO  [mixer_module] allocateControl time: 44 
INFO  [mixer_module] DPscaled_LPCA time: 68 
INFO  [mixer_module] DP_LPCA time: 188 
INFO  [mixer_module] allocator_dir_LPwrap_4 time: 26 
## for compare, the running time of older version dir_alloc_sim which generated by matlab and general inversion method:
INFO  [mixer_module] dir_alloc_sim time: 438 
INFO  [mixer_module] inv time: 2 
# for attainable m_higher {0.0,  0.0,  10.0f}
INFO  [mixer_module] DP_LPCA_prio time: 158
# for unattainable m_higher {0.0,  0.0,  60.0f}, have to restoring
INFO  [mixer_module] DP_LPCA_prio time: 250 
```
output for attainable m_higher {0.0,  0.0,  0.0f} is:
```sh
INFO  [mixer_module] allocateControl: u1: -0.109583, u2: -0.234828, u3: 0.349100, u4: -0.004688. 

INFO  [mixer_module] allocateControl time: 52 

INFO  [mixer_module] DPscaled_LPCA: u1: -0.229341, u2: -0.115069, u3: 0.229341, u4: 0.115069. 

INFO  [mixer_module] DPscaled_LPCA time: 96 

INFO  [mixer_module] DP_LPCA: u1: -0.229341, u2: -0.115069, u3: 0.229341, u4: 0.115069. 

INFO  [mixer_module] DP_LPCA time: 227 

INFO  [mixer_module] DP_LPCA_prio: u1: -0.229341, u2: -0.115069, u3: 0.229341, u4: 0.115069. 

INFO  [mixer_module] DP_LPCA_prio time: 240
```

with restoring, add  40us or 50us.

output for attainable m_higher {0.0,  0.0,  10.0f} is:
```sh
INFO  [mixer_module] allocateControl: u1: 0.044158, u2: 0.196098, u3: 0.349100, u4: 0.349100. 

INFO  [mixer_module] allocateControl time: 55 

INFO  [mixer_module] DPscaled_LPCA: u1: 0.044158, u2: 0.196098, u3: 0.349100, u4: 0.349100. 

INFO  [mixer_module] DPscaled_LPCA time: 89 

INFO  [mixer_module] DP_LPCA: u1: 0.044158, u2: 0.196098, u3: 0.349100, u4: 0.349100. 

INFO  [mixer_module] DP_LPCA time: 232 

INFO  [mixer_module] DP_LPCA_prio: u1: 0.349100, u2: 0.349100, u3: 0.349100, u4: 0.349100. 

INFO  [mixer_module] DP_LPCA_prio time: 354
```

ToDo: delete inital of A b h