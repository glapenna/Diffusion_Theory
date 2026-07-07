#!/bin/bash -f

BINDIR=${HOME}/diffusion/bin

module purge
module load profile/base
module load openblas/0.3.24--gcc--12.2.0
module load openmpi/4.1.6--gcc--12.2.0                           
module load netlib-scalapack/2.2.0--openmpi--4.1.6--gcc--12.2.0  

gfortran -c solveSF.f90
gfortran -o ${BINDIR}/solveSF solveSF.o -L${NETLIB_SCALAPACK_LIB}  -lopenblas -lscalapack
