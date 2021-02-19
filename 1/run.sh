#!/bin/bash
echo "Running process with $1 $2"
MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin
$MPI_PATH/mpicxx -Wunused-variable main$1.cpp -o main$1 -std=c++11
$MPI_PATH/mpirun -machinefile $PBS_NODEFILE -np $2 ./main$1