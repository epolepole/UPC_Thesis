#!/bin/bash
mkdir build
mkdir build/bin
cp SetUpData.ini build/bin/
echo "Coying Mesh_$3 to build/bin/Mesh"
cp -r Mesh_$3/ build/bin/Mesh
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
    build_type=Debug
    save_iter=0
    echo Running without iterations

elif [ $1 == "P" ] ; then
    verbose_makefile=OFF
    build_type=Release
    save_iter=0
    print_info=1
    echo Running with prints

elif [ $1 == "I" ] ; then
    verbose_makefile=OFF
    build_type=Release
    save_iter=1
    echo Runing with iterations

elif [ $1 == "T" ] ; then
    verbose_makefile=ON
    build_type=Debug
    save_iter=1
    echo Tracing

elif [ $1 == "D" ] ; then
    verbose_makefile=ON
    build_type=Debug
    save_iter=1
    print_info=1
    echo Debugging

elif [ $1 == "GDB" ] ; then
    verbose_makefile=ON
    build_type=Debug
    save_iter=1
    print_info=1
    echo Debugging with GDB
fi

cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=$verbose_makefile -DCMAKE_BUILD_TYPE=$build_type\
 -DSAVE_ITER=$save_iter -DNN=$3 -DNM=$3 -DNL=$3 -DLAT=$4 -DPRINT_INFO=$print_info
make
cd ../
