#!/bin/bash
clear
./compile.sh

t_end=2.5001
uemax=0.5

mul=5.00e-3
mug=2.50
ul=0.5
ug=53.33333
theta0=70

#mpirun ./blur3D-ibm-mpi 1> out 2> log
#mpiexec -np 20 ./blur-ibm-openmp $t_end $mul $mug $ul $ug $theta0 $uemax 1> out 2> log

./blur-ibm-openmp $t_end $mul $mug $ul $ug $theta0 $uemax 1> out 2> log
