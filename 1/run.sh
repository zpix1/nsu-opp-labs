#!/bin/bash
mpic++ main.cpp -o main
mpirun -np 3 ./main