#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <fstream>

#include <mpi.h>


#define DEBUG(var) \
            do { std::cout << p_rank << " has " << #var << ": " << var << std::endl; } while (0)

void fill(double* x, int N, double value) {
    for (int i = 0; i < N; i++) {
        x[i] = value;
    }
}

void mulAB(int N1, int N2, int N3, double* A, double* B, double* C)  {
    // ..
}

void split(int* starts, int* sizes, int N, int pc, int k) {
    int offset = 0;
    for (int i = 0; i < pc; i++) {
        sizes[i] = N / pc;
        if (i < N % pc) {
            sizes[i]++;
        }

        sizes[i] *= k;

        starts[i] = offset;
        offset += sizes[i];
    }
}

void mpi_mat_mat_mul(int m, int n, int k, double* A, double* B, double* C, MPI_Comm comm2d, int* dims, int* periods) {
    int coords[2];
    MPI_Cart_get(comm2d, 2, dims, periods, coords);
    int ranky = coords[0];
    int rankx = coords[1];

    if (ranky == 0) {
        // Scatter
    }
    if (rankx == 0) {
        // Scatter
    }

    // BCast

}

int main(int argc, char** argv) {    
    int p_rank;
    int p_count;

    const int N1 = 5;
    const int N2 = 6;
    const int N3 = 7;
    const int P1 = 2;
    const int P2 = 2;

    double* matrix_A;
    double* matrix_B;
    double* matrix_C;
    
    double start, end;

    // INIT MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);

    assert(P1 * P2 == p_count);

    // LOAD MATRIX
    if (p_rank == 0) {
        matrix_A = new double[N1*N2];
        fill(matrix_A, N1*N2, 5.);
        matrix_B = new double[N2*N3];
        fill(matrix_B, N2*N3, 7.);
        matrix_C = new double[N2*N2];
    }

    int dims[2] = {0, 0}, periods[2] = {0, 0};
    MPI_Comm comm2d;
    MPI_Dims_create(p_count, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm2d);

    mpi_mat_mat_mul(N1, N3, N2, matrix_A, matrix_B, matrix_C, comm2d, dims, periods);

}