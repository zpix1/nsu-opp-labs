#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <fstream>

#include <mpi.h>
#include <mpe.h>

int p_rank;
int p_count;

#define DEBUG(var) \
            do { std::cout << p_rank << " has " << #var << ": " << var << std::endl; } while (0)

const double EPS = 1e-10;

double scalar(const double* a, const double* b, int N) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += a[i] * b[i];
    }
    return res;
}

double vscalar(const double* a, const double* b, int N, int* v_starts, int* v_sizes, int p_rank) {
    double res = 0;
    for (int i = 0; i < v_sizes[p_rank]; i++) {
        res += a[v_starts[p_rank] + i] * b[v_starts[p_rank] + i];
    }
    double ans;
    MPI_Allreduce(&res, &ans, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return ans;
}

void dump(double* x, int N, std::string fname) {
    std::ofstream f(fname);
    for (int i = 0; i < N; i++) {
        f << x[i] << " ";
    }
    f << std::endl;
    f.close();
}

int main(int argc, char** argv) {
    double start, end;
    int N;

    int evtid_begin_matrix_mul, evtid_end_matrix_mul;
    int evtid_begin_addition_calculations, evtid_end_addition_calculations;
    int evtid_begin_allgv, evtid_end_allgv;
    // === Initialize MPI ===

    MPI_Init(&argc, &argv);
    MPE_Init_log();

    evtid_begin_matrix_mul = MPE_Log_get_event_number();
    evtid_end_matrix_mul = MPE_Log_get_event_number();
    evtid_begin_addition_calculations = MPE_Log_get_event_number();
    evtid_end_addition_calculations = MPE_Log_get_event_number();
    evtid_begin_allgv = MPE_Log_get_event_number();
    evtid_end_allgv = MPE_Log_get_event_number();

    MPE_Describe_state(evtid_begin_matrix_mul, evtid_end_matrix_mul, "Vector-Matrix Mul", "red");
    MPE_Describe_state(evtid_begin_addition_calculations, evtid_end_addition_calculations, "Other vector calculations", "blue");
    MPE_Describe_state(evtid_begin_allgv, evtid_end_allgv, "MPI_Allgatherv", "yellow");

    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    // === Load data and split data across processes ===
    double* A;
    double* b;
    double* x;

    if (p_rank == 0) {
        std::ifstream f("input.dat");
        f >> N;
        A = new double[N*N];
        b = new double[N];
        x = new double[N*N];
        for (int i = 0; i < N*N; i++) {
            f >> A[i];
        }
        for (int i = 0; i < N; i++) {
            f >> b[i];
            x[i] = 0;
        }
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    const int MAX_PART_SIZE = N / p_count + 1;

    double* A_part = new double[MAX_PART_SIZE * N];
    double* b_part = new double[MAX_PART_SIZE];
    double* x_part = new double[MAX_PART_SIZE];

    int* b_starts = new int[p_count];
    int* b_sizes = new int[p_count];
    int* A_starts = new int[p_count];
    int* A_sizes = new int[p_count];

    int b_offset = 0;
    for (int i = 0; i < p_count; i++) {
        b_sizes[i] = N / p_count;
        // To remove any gaps
        if (i < N % p_count) {
            b_sizes[i]++;
        }
        b_starts[i] = b_offset;
        b_offset += b_sizes[i];
        // By rows
        A_starts[i] = b_starts[i] * N;
        A_sizes[i] = b_sizes[i] * N;
    }

    std::cout << "process " << p_rank << " will do " << b_starts[p_rank] << "-" << b_starts[p_rank]+b_sizes[p_rank]  << std::endl;

    MPI_Scatterv(A, A_sizes, A_starts, MPI_DOUBLE, A_part, A_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b, b_sizes, b_starts, MPI_DOUBLE, b_part, b_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(x, b_sizes, b_starts, MPI_DOUBLE, x_part, b_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* r = new double[N];
    double* z = new double[N];
    double* Az = new double[N];
    double* part_buf0 = new double[b_sizes[p_rank]];
    double* part_buf1 = new double[b_sizes[p_rank]];
    double* full_buf = new double[N];

    if (p_rank == 0) {
        start = MPI_Wtime();
    }


    // 1 - initialization
    // r0 = b0 0- Ax0
    // z0 = r0
    double r_square = 0;
    double b_square = 0;
    MPE_Log_event(evtid_begin_matrix_mul, p_rank, (char*)0);
    double b_square_part = 0;
    for (int i = 0; i < b_sizes[p_rank]; i++) {
        part_buf0[i] = 0;
        for (int j = 0; j < N; j++) {
            part_buf0[i] += A_part[i * N + j] * x_part[j];
        }
        part_buf0[i] = b_part[i] - part_buf0[i];
        b_square_part += b_part[i] * b_part[i];
    }
    MPE_Log_event(evtid_end_matrix_mul, p_rank, (char*)0);

    MPI_Allreduce(&b_square_part, &b_square, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allgatherv(part_buf0, b_sizes[p_rank], MPI_DOUBLE, r, b_sizes, b_starts, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(part_buf0, b_sizes[p_rank], MPI_DOUBLE, z, b_sizes, b_starts, MPI_DOUBLE, MPI_COMM_WORLD);
    
    r_square = vscalar(r, r, N, b_starts, b_sizes, p_rank);

    // 2 - main cycle
    // int i = 0;
    int converges = false;
    while (!converges) {
        MPE_Log_event(evtid_begin_matrix_mul, p_rank, (char*)0);
        for (int i = 0; i < b_sizes[p_rank]; i++) {
            part_buf0[i] = 0;
            for (int j = 0; j < N; j++) {
                part_buf0[i] += A_part[i * N + j] * z[j];
            }
        }
        MPE_Log_event(evtid_end_matrix_mul, p_rank, (char*)0);
        MPE_Log_event(evtid_begin_allgv, p_rank, (char*)0);
        MPI_Allgatherv(part_buf0, b_sizes[p_rank], MPI_DOUBLE, Az, b_sizes, b_starts, MPI_DOUBLE, MPI_COMM_WORLD);
        MPE_Log_event(evtid_end_allgv, p_rank, (char*)0);
        MPE_Log_event(evtid_begin_addition_calculations, p_rank, (char*)0);

        double alpha = r_square / vscalar(Az, z, N, b_starts, b_sizes, p_rank);
        
        for (int i = 0; i < b_sizes[p_rank]; i++) {
            x_part[i] = x_part[i] + alpha * z[b_starts[p_rank] + i];
        }

        // full_buf - z_n+1
        for (int i = 0; i < N; i++) {
            full_buf[i] = r[i] - alpha * Az[i];
        }
        double r_new_square = vscalar(full_buf, full_buf, N, b_starts, b_sizes, p_rank);
        double beta = r_new_square / r_square;
        for (int i = 0; i < N; i++) {
            r[i] = full_buf[i];
        }
        
        for (int i = 0; i < N; i++) {
            z[i] = r[i] + beta * z[i];
        }
        if (sqrt(r_square/b_square) < EPS) {
            converges = true;
        }
        r_square = r_new_square;
        MPE_Log_event(evtid_end_addition_calculations, p_rank, (char*)0);
    }
    MPI_Gatherv(x_part, b_sizes[p_rank], MPI_DOUBLE, x, b_sizes, b_starts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPE_Finish_log("main2mpe.clog2");
    if (p_rank == 0) {
        end = MPI_Wtime();
        std::cout << "DONE: " << p_count << ": " << end - start << std::endl;
        dump(x, N, "output.dat");
    }

    if (p_rank == 0) {
        delete[] A;
        delete[] b;

        delete[] x;
    }

    MPI_Finalize();

    delete[] z;
    delete[] Az;
    delete[] part_buf0;
    delete[] part_buf1;
    delete[] full_buf;
    delete[] b_starts;
    delete[] b_sizes;
    delete[] A_starts;
    delete[] A_sizes;

    return 0;
}