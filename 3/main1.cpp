#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <fstream>

#include <mpi.h>

int p_rank;
int p_count;

#define DEBUG(var) \
            do { std::cout << p_rank << " has " << #var << ": " << var << std::endl; } while (0)

void fill(double* x, int N, double value) {
    for (int i = 0; i < N; i++) {
        x[i] = value;
    }
}

int main(int argc, char** argv) {
    double start, end;
    
    // === Initialize MPI ===

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    
    // === Load data and split data across processes ===
    
    MPI_Finalize();
    return 0;
}