#!/bin/bash
mpic++ main.cpp -o main
mpirun -np 2 ./main