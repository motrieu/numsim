#!/bin/bash

module use /usr/local.nfs/sgs/modulefiles
module load gcc/10.2
module load openmpi/3.1.6-gcc-10.2
module load vtk/9.0.1
module load cmake/3.18.2

cmake ../..
cd ..
make install
cd src

srun -n 16 ./numsim lid_driven_cavity.txt
