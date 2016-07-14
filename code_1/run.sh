#!/bin/bash

cd build/bin
if [ $1 == "R" ] ; then
    upcrun -shared-heap 256MB -backtrace LBMSolver
else
    if [ $1 == "I" ] ; then
        upcrun -shared-heap 256MB -backtrace LBMSolver
    else
        upcrun -freeze=0 -backtrace LBMSolver
    fi
fi
cd ../../
