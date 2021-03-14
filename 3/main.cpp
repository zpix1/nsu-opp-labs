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

void mpi_mat_mat_mul(int m, int n, int k, double* A, double* B, double* C, MPI_Comm comm, int* p) {
    int size;
    MPI_Comm_size(comm, &size);

    int periods[2] = {0, 0};
    MPI_Comm comm2d;
    MPI_Dims_create(size, 2, p);
    MPI_Cart_create(comm, 2, p, periods, 0, &comm2d);

    int rank;
    MPI_Comm_rank(comm2d, &rank);

    int coords[2];
    MPI_Cart_get(comm2d, 2, p, periods, coords);

    // ?
    MPI_Comm comm1d[2];
    int remains[2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            remains[j] = (i == j);
            MPI_Cart_sub(comm2d, p, &comm1d[i]);
        }
    }

    int nn[2];

    // should be divisible without remainder
    nn[0] = m / p[0];
    nn[1] = n / p[1];

    double* AA = new double[nn[0] * k];
    double* BB = new double[nn[1] * k];
    double* CC = new double[nn[0] * nn[1]];

    if (coords[1] == 0) {
        MPI_Scatter(A, k * nn[0], MPI_DOUBLE, AA, k * nn[0], MPI_DOUBLE, 0, comm1d[0]);
    }

    MPI_Datatype vector_t;
    MPI_Datatype resized_vector_t;
    MPI_Datatype recv_t;
    
    int *dispc, *countc;
    if (coords[0] == 0) {
        MPI_Type_vector(k, nn[1], n, MPI_DOUBLE, &vector_t);
        MPI_Type_commit(&vector_t);

        MPI_Type_create_resized(vector_t, 0, nn[1] * sizeof(double), &resized_vector_t);
        MPI_Type_commit(&resized_vector_t);

        MPI_Type_contiguous(k * nn[1], MPI_DOUBLE, &recv_t);
        MPI_Type_commit(&recv_t);

        MPI_Scatter(B, 1, resized_vector_t, BB, 1, recv_t, 0, comm1d[1]);

        dispc = new int[p[0] * p[1]];
        countc = new int[p[0] * p[1]];

        for (int i = 0; i < p[0]; i++) {
            for (int j = 0; j < p[1]; j++) {
                dispc[i*p[1] + j] = (i * p[1] * nn[0] + j);
                countc[i*p[1] + j] = 1;
            }
        }
    }
    
    MPI_Bcast(AA, k * nn[0], MPI_DOUBLE, 0, comm1d[1]);
    MPI_Bcast(BB, k * nn[1], MPI_DOUBLE, 0, comm1d[0]);

    #define AA(i,j) AA[k*i + j]
    #define BB(i,j) BB[nn[1]*i + j]
    #define CC(i,j) CC[nn[1]*i + j]
    
    for (int i = 0; i < nn[0]; i++) {
        for (int j = 0; j < nn[1]; j++) {
            CC(i, j) = 0.0;
            for (int g = 0; g < k; g++) {
                CC(i, j) += AA(i, g) * BB(g, j);
            }
        }
    }

    // if (rank == 3) {
    //     for (int i = 0; i < nn[0] * nn[1]; i++) {
    //         printf("%4f ", CC[i]);
    //     }
    //     printf("\n");
    // }

    MPI_Datatype recv_vector_t, send_vector_t, resized_recv_vector_t;
    
    MPI_Type_contiguous(nn[0] * nn[1], MPI_DOUBLE, &send_vector_t);
    MPI_Type_commit(&send_vector_t);

    MPI_Type_vector(nn[0], nn[1], n, MPI_DOUBLE, &recv_vector_t);
    MPI_Type_commit(&recv_vector_t);

    MPI_Type_create_resized(recv_vector_t, 0, nn[1] * sizeof(double), &resized_recv_vector_t);
    MPI_Type_commit(&resized_recv_vector_t);

    MPI_Gatherv(CC, 1, send_vector_t, C, countc, dispc, resized_recv_vector_t, 0, comm2d);
}

int main(int argc, char** argv) {    
    int p_rank;
    int p_count;

    const int N1 = 8;
    const int N2 = 8;
    const int N3 = 8;
    const int P1 = 1;
    const int P2 = 1;

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
        fill(matrix_A, N1*N2, 0.);
        for (int i = 0; i < N1; i++) {
                matrix_A[i * N2 + i] = 2.0;
        }
        matrix_B = new double[N2*N3];
        fill(matrix_B, N2*N3, 7.);
        matrix_C = new double[N2*N2];
        printf("Matrix loading done;\n");
    }

    int p[2] = {P1, P2};

    mpi_mat_mat_mul(N1, N3, N2, matrix_A, matrix_B, matrix_C, MPI_COMM_WORLD, p);

    if (p_rank == 0) {
        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N3; j++) {
                printf("%4.0f ", matrix_C[N3*i + j]);
            }
            printf("\n");
        }
    }

    MPI_Finalize();

    return 0;
}