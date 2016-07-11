#!/bin/bash
mkdir build
mkdir build/bin
cp SetUpData.ini build/bin/
cp -r Mesh/ build/bin/
cd build
#rm -fr *

if [ -z "$2" ] ; then
threads=8
else
threads=$2
fi


if [ $1 == "R" ] ; then
    echo Running without iterations
    cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release
    make
else
    if [ $1 == "I" ] ; then
        echo Runing with iterations
        cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release -DSAVE_ITER=1
        make
    else
        echo Debugging
        cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Debug -DSAVE_ITER=1
        make verbose=1
    fi
fi

cd ../
#cmake .. -DCMAKE_VERBOSE_MAKEFILE=OFF