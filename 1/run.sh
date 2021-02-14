#!/bin/bash
mpic++ -Wunused-variable main$1.cpp -o main
mpirun -np $2 ./main