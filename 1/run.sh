#!/bin/bash
mpic++ -Wunused-variable main.cpp -o main
mpirun -np 3 ./main