#UPC_Thesis

README FILE

The code is on the folder code_1 and the scripts to run it can be found there.


This code already includes some meshes to test with.


In order to compile the code, the berkeley upc runtime must be installed, and the bin and include folders, must be aded to the PATH variable.



There are two options for compiling and running the code. With CMake, or withouth CMake.

Both uses their correspondent script (go.sh or go_noCMake.sh)

The imput parameters for both scripts are:
1: R/I/P/T/D/GDB
	R: release
	I: saving iterations every certain number
	P: printing extra info on screen
	T: run with backtrace debugging
	D: run with debbugging
	GDB: run with gdb and freezing the code for linking the gdb to the process

2: Number of threads (In the case of the cubes, the threads must be a power of 3 (1,8,27,64...))
3: Size of the mesh. Select the size from the available meshes on the folder, and taking into account the number of threads. The size must be divisible by the cubic root of the number of threads. (With 27 threads, the size must be multiple of 3)
4: Result of this previous division (Threads=27, Size=63 -> 21)


Examples:

./go.sh R 8 16 8
./go_noCMake.sh I 27 63 21





/////
CMAKE

Important, before using CMake, an internal CMake option must be deactivated in order to use the upc compiler. The file to edit is located on the folder:
	/usr/share/cmake-3.5/Modules/Platform/

The file to edit is the one refering to your C compiler, for example, if using the Linux GNU compiler, the file to edit is:
	Linux-GNU.cmake

On this file, one has to move to the line that sais:
	set(CMAKE_SHARED_LIBRARY_LINK_\$\{lang\}_FLAGS "-rdynamic")

and change it for:
	set(CMAKE_SHARED_LIBRARY_LINK_\$\{lang\}_FLAGS "")

once that is done, the code can be compiled and run at the same time with the script
	go.sh



////////
NO CMake
The code can be simply compiled with the script named:
	go_noCMake.sh







The version here is just the cubes approach, the repository with all the changes and branches, with the layers approach can be found at:
	git@github.com:epolepole/UPC_Thesis.git