#!/bin/bash
echo "Running from"
pwd
cd build/bin
if [ $1 == "R" ] ; then
    upcrun -shared-heap 1000MB LBMSolver
elif [ $1 == "P" ] ; then
    upcrun -shared-heap 1000MB LBMSolver
elif [ $1 == "I" ] ; then
    upcrun -shared-heap 1000MB LBMSolver
elif [ $1 == "T" ] ; then
    upcrun -shared-heap 1000MB -trace LBMSolver
elif [ $1 == "D" ] ; then
    upcrun -shared-heap 1000MB -backtrace  LBMSolver
elif [ $1 == "GDB" ] ; then
    upcrun -shared-heap 1000MB -freeze=0 -backtrace LBMSolver

fi
cd ../../
