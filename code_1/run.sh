#!/bin/bash

cd build/bin
if [ $1 == "R" ] ; then
    upcrun -backtrace LBMSolver
else
    upcrun -freeze=0 -backtrace LBMSolver
fi
cd ../../