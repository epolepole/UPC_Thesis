#!/bin/bash

#1: R I D (Release, Iterations, Debug)
#2: Number of threads
#3: Mesh size

#./clean.sh
./rebuild.sh $1 $2 $3 $4
./run.sh $1