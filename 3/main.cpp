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

void split_matrix_parts(int* starts, int* sizes, int N, int pc) {
    int offset = 0;
    for (int i = 0; i < pc; i++) {
        sizes[i] = N / pc;
        if (i < N % pc) {
            sizes[i]++;
        }
        
        sizes[i] *= N;

        starts[i] = offset;
        offset += sizes[i];
    }
}

int main(int argc, char** argv) {
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

    // SPLIT DATA
    int* p1_starts = new int[P1];
    int* p1_sizes = new int[P1];

    int* p2_starts = new int[P2];
    int* p2_sizes = new int[P2];

    int p1_offset = 0;
    int p2_offset = 0;

    split_matrix_parts(p1_starts, p1_sizes, N2, P1);
    split_matrix_parts(p2_starts, p2_sizes, N2, P2);

    // SCATTER DATA
    double* matrix_A_part = new double[p1_sizes[p_rank] * N2];
    double* matrix_B_part = new double[p2_sizes[p_rank] * N2];

    MPI_Scatterv(matrix_A, p1_sizes, p1_starts, MPI_DOUBLE, matrix_A_part, p1_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(matrix_B, p2_sizes, p2_starts, MPI_DOUBLE, matrix_B_part, p2_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // MULTIPLY
    
}