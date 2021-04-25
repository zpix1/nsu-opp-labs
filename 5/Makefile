all:
	mpic++ -pthread -std=c++14 -g -Wall main.cpp -o main
run:
	mpirun -np 4 ./main
mpe:
	mpecc -mpicc='/opt/intel/impi/5.0.1.035/intel64/bin/mpicxx' -std=c++11 -mpilog -o mainmpe main.cpp
mperun:
	mpirun -trace -np 4 ./mainmpe