#!/bin/bash
mkdir build
mkdir build/bin
cp SetUpData.ini build/bin/
cp -r Mesh/ build/bin/
cd build
#rm -fr *
cmake .. -DUPC=1 -DTHREADS=8 -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release
#cmake .. -DCMAKE_VERBOSE_MAKEFILE=OFF
make verbose=1