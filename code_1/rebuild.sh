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
    cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release
    make
else
    cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Debug
    make verbose=1
fi

cd ../
#cmake .. -DCMAKE_VERBOSE_MAKEFILE=OFF