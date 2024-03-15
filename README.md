## About
This is a repos. for test new control allocation algorithm developed from [ simulation of the book "aircraft control allocation"](https://github.com/mengchaoheng/aircraft-control-allocation) and [qcat](https://github.com/mengchaoheng/qcat).

The new algorithm will used by a ductedfan UAV.

![image](https://user-images.githubusercontent.com/43166007/143251407-38f54ce4-739e-40d5-85ff-8961ab41fff8.png)

![image](https://user-images.githubusercontent.com/43166007/143251429-d33c9e1a-9cb0-4889-9e8b-cabe5bf9184a.png)

1. the lib. of qcat and aircraft-control-allocation-book-simulation is used as submodules here.

2. the ac.slx and some files from the compare of allocation in Plan-D (`alloc_compa`).


3. test.m, constraint.m and twin.slx come from `allocation m2C`，view the attainable set。

4. some new algorithm is show here.

# Some files and folders
submodules used sa lib:
* aircraft-control-allocation-book-simulation
* qcat

new algorithm
* LP
* PCA
* QP

simulink test and some test:
* `ac.slx` and `twin.slx`
* test_xxx.m

come from the old version modified `aircraft-control-allocation-book-simulation`:
* some_modified_function 
 
come from the old version `control allocation`:
* s-function_used_in PlanD 
* plot_actuator.m 


