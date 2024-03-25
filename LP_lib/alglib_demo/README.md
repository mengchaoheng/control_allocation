# How to build

## linux
For the demo example
### using commander line tools

run following:
1. create and cd to build folder:
```
mkdir build && cd build
```
2. build the code by g++:
```
$ g++ -I ../src/alglib-cpp/src/ -o demo.out ../src/demo/*.cpp ../src/alglib-cpp/src/*.cpp
```
3. run 
```
./demo.out
```

If create new target, change the name  'demo'. For dir:
```Console
g++ -I ../src/alglib-cpp/src/ -o minlp_basic.out ../src/allocator_dir/*.cpp  ../src/alglib-cpp/src/*.cpp -w -O3
```
then run:
```Console
./minlp_basic.out
```


## test
demo
```Console
➜  build git:(main) ✗ g++ -I ../src/alglib-cpp/src/ -o demo.out ../src/demo/*.cpp ../src/alglib-cpp/src/*.cpp -w
➜  build git:(main) ✗ ./demo.out 
Performance is 2.1 GFLOPS
```
-O3
```Console
➜  build git:(main) ✗ g++ -I ../src/alglib-cpp/src/ -o demo.out ../src/demo/*.cpp ../src/alglib-cpp/src/*.cpp -w -O3
➜  build git:(main) ✗ ./demo.out 
Performance is 14.2 GFLOPS
```
test_c
```Console
➜  build git:(main) g++ -I ../src/alglib-cpp/src/ -o test_c.out  ../src/alglib-cpp/src/*.cpp ../src/alglib-cpp/tests/test_c.cpp  -w 

➜  build git:(main) ✗ ./test_c.out
SEED: 1711070688
COMPILER: GCC
HARDWARE: 64-bit
HARDWARE: little-endian
CPU:   unknown
CORES: 1 (serial version)
LIBS:  
CPUID:
OS: unknown
TESTING MODE: single core
ablasf                           OK
hqrnd                            OK
a
```
test_LPWrap, 
```Console
$ g++ -I ../src/alglib-cpp/src/ -o LPWrap_test.out  ../src/LPWrap/*.cpp ../src/LPWrap/*.c ../src/alglib-cpp/src/*.cpp -w  -O3

$ ./LPWrap_test.out
```



## other
为了使用g++编译C++项目并将生成的文件保存到新的build文件夹，你可以在项目根目录下创建一个构建脚本或直接在命令行中执行以下步骤：

创建build文件夹：

mkdir build

进入build文件夹：

cd build

使用g++编译源代码，例如：

g++ -std=c++11 -o my_program main.cpp

这里假设你的源文件是main.cpp，生成的可执行文件名为my_program。-std=c++11是用来指定C++标准的，你可以根据需要替换为其他标准，如c++14, c++17, c++20等。

如果你有多个源文件，你可以将它们全部列出，或者使用通配符：

g++ -std=c++11 -o my_program *.cpp

如果你想要自动化这个过程，你可以写一个简单的构建脚本，例如build.sh（在Unix-like系统中）：

#!/bin/bash
 
# 清理旧的build文件夹
rm -rf build
mkdir build
 
# 进入build文件夹并编译
cd build && g++ -std=c++11 -o my_program ../main.cpp && cd ..

确保给这个脚本可执行权限：

chmod +x build.sh

然后运行它：

./build.sh

这个脚本会清理旧的build文件夹，创建一个新的，并编译源代码文件到build文件夹中。


##  C++和C混合编译
在C++中混合编译C代码时，需要使用extern "C"来告知编译器这部分代码应当使用C的链接方式，而不是C++的。这是因为C和C++对函数名字的修饰(name mangling)方式不同。

例如，假设你有一个C函数int c_function()在C文件中定义，你想在C++文件中调用它。

C文件 c_code.c:
```c
#include <stdio.h>
 
int c_function() {
    printf("Called C function\n");
    return 0;
}
```
C++文件 cpp_code.cpp:
```Cpp
extern "C" {
    int c_function(); // 声明C函数，注意不要包含分号
}
 
int main() {
    c_function(); // 调用C函数
    return 0;
}
```
在这个例子中，extern "C"告诉编译器c_function是一个在C语言链接命名规则下编译的函数。在编译和链接时，C++代码会正确地链接到C函数。

