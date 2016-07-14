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

if [ $1 == "R" ] ; then
    echo Running without iterations
    cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release -DNN=$3 -DNM=$3 -DNL=$3
    make
else
    if [ $1 == "I" ] ; then
        echo Runing with iterations
        cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_BUILD_TYPE=Release -DSAVE_ITER=1 -DNN=$3 -DNM=$3 -DNL=$3
        make
    else
        echo Debugging
        cmake .. -DUPC=1 -DTHREADS=$threads -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Debug -DSAVE_ITER=1 -DNN=$3 -DNM=$3 -DNL=$3
        make verbose=1
    fi
fi

cd ../