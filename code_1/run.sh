#!/bin/bash
echo "Running from"
pwd
cd build/bin
if [ $1 == "R" ] ; then
    upcrun -shared-heap 1024MB -backtrace LBMSolver
elif [ $1 == "P" ] ; then
    upcrun -shared-heap 1024MB LBMSolver
elif [ $1 == "I" ] ; then
    upcrun -shared-heap 1024MB LBMSolver
elif [ $1 == "T" ] ; then
    upcrun -shared-heap 1024MB -trace LBMSolver
elif [ $1 == "D" ] ; then
    upcrun -shared-heap 1024MB -backtrace  LBMSolver
elif [ $1 == "GDB" ] ; then
    upcrun -shared-heap 1024MB -freeze=0 -backtrace LBMSolver

fi
cd ../../
