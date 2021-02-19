#!/bin/bash

#PBS -l select=2:ncpus=8:mpiprocs=16:mem=2000m,place=free
#PBS -l walltime=01:00:00
#PBS -q S3077545
# cp runner.py settings.py generate_input.py main*.cpp run.sh $PBS_O_WORKDIR

cd $PBS_O_WORKDIR
echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin
$MPI_PATH/mpicxx -Wunused-variable main1.cpp -o main1 -std=c++11
$MPI_PATH/mpirun -machinefile $PBS_NODEFILE -np 1 ./main1