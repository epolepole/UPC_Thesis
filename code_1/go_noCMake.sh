#!/bin/bash
mkdir build
mkdir build/bin
cp -f SetUpData.ini build/bin/
echo "Coying Mesh_$3 to build/bin/Mesh"
cp -fr Mesh_$3/ build/bin/Mesh
cd build
#rm -fr *

if [ -z "$2" ] ; then
threads=8
else
threads=$2
fi

rm -fr bin/Results

print_info=0


if [ $1 == "R" ] ; then
    verbose_makefile=OFF
    build_type=Release
    gdb=""
    save_iter=0
    echo Running without iterations

elif [ $1 == "P" ] ; then
    verbose_makefile=OFF
    build_type=Release
    gdb=""
    save_iter=0
    print_info=1
    echo Running with prints

elif [ $1 == "I" ] ; then
    verbose_makefile=OFF
    build_type=Release
    gdb=""
    save_iter=1
    echo Runing with iterations

elif [ $1 == "T" ] ; then
    verbose_makefile=ON
    build_type=Debug
    gdb="-g"
    save_iter=1
    echo Tracing

elif [ $1 == "D" ] ; then
    verbose_makefile=ON
    build_type=Debug
    gdb="-g"
    save_iter=1
    print_info=1
    echo Debugging

elif [ $1 == "GDB" ] ; then
    verbose_makefile=ON
    build_type=Debug
    gdb="-g"
    save_iter=1
    print_info=1
    echo Debugging with GDB


fi


flags="-DLAT=$4 -DNL=$3 -DNM=$3 -DNN=$3 -I../include -T $2  -D__SAVE_ITER__=$save_iter -D__PRINT_INFO__=$print_info"
upcc  $flags $gdb -o main.c.o                    -c ../src/main.c
upcc  $flags $gdb -o BoundaryConditions.c.o      -c ../src/BoundaryConditions.c
upcc  $flags $gdb -o CellFunctions.c.o           -c ../src/CellFunctions.c
upcc  $flags $gdb -o ComputeResiduals.c.o                    -c ../src/ComputeResiduals.c
upcc  $flags $gdb -o FilesReading.c.o                    -c ../src/FilesReading.c
upcc  $flags $gdb -o FilesWriting.c.o                    -c ../src/FilesWriting.c
upcc  $flags $gdb -o Iterate.c.o                    -c ../src/Iterate.c
upcc  $flags $gdb -o tests.c.o                    -c ../src/tests.c
upcc  $flags $gdb -o ShellFunctions.c.o                    -c ../src/ShellFunctions.c

upcc -T $2  -D__SAVE_ITER__=$save_iter  main.c.o ShellFunctions.c.o ComputeResiduals.c.o tests.c.o FilesWriting.c.o FilesReading.c.o CellFunctions.c.o BoundaryConditions.c.o Iterate.c.o  -o bin/LBMSolver  -lm 

cd ../




./run.sh $1