#!/bin/bash

#PBS -l select=2:ncpus=8:mpiprocs=16:mem=2000m,place=free
#PBS -l walltime=01:01:00

cd $PBS_O_WORKDIR

mpirun -machinefile $PBS_NODEFILE -np 16 ./mainmpe 2000 2000 2000 4 4