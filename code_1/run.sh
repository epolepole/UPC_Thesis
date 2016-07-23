#!/bin/bash
echo "Running from"
pwd
cd build/bin
if [ $1 == "R" ] ; then
    upcrun -shared-heap 512MB LBMSolver
elif [ $1 == "I" ] ; then
    upcrun -shared-heap 512MB LBMSolver
elif [ $1 == "T" ] ; then
    upcrun -shared-heap 512MB -trace LBMSolver
elif [ $1 == "D" ] ; then
    upcrun -freeze=0 -backtrace -backtrace -trace LBMSolver

fi
cd ../../
