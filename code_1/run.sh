#!/bin/bash

cd build/bin
if true ; then
upcrun -backtrace LBMSolver
else
upcrun -freeze=0 -backtrace LBMSolver
fi