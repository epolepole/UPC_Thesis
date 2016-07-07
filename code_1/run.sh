#!/bin/bash

cd build/bin
if [ $1 == "R" ] ; then
    upcrun -backtrace LBMSolver
else
    upcrun  -backtrace LBMSolver
fi
cd ../../