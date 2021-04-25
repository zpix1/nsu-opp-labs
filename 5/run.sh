#!/bin/bash

# NPC=mpiprocs=select*ncpus
# 1 <= select <= 2

#PBS -l select=1:ncpus=4:mpiprocs=4:mem=4000m,place=free:exclhost
#PBS -l walltime=01:01:00

export MPE_HOME=$HOME/mpe2-nusc-built
export PATH=$PATH:$MPE_HOME/bin

cd $PBS_O_WORKDIR

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin

mpecc -mpicc='/opt/intel/impi/5.0.1.035/intel64/bin/mpicxx' -std=c++11 -mpilog -o mainmpe main.cpp
mpirun -machinefile $PBS_NODEFILE -trace -np 4 ./mainmpe

# $MPI_PATH/mpicxx main.cpp -o emain1 -std=c++11 -O2 -xhost=native
# $MPI_PATH/mpirun -machinefile $PBS_NODEFILE -np 1 ./emain1
