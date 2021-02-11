#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>

#include <mpi.h>

const double EPS = 1e-5;

void mulMV(const double* A, const double* x, int N, double* res)  {
    for (int i = 0; i < N; i++) {
        res[i] = 0;
        for (int j = 0; j < N; j++) {
            res[i] += A[i * N + j] * x[j];
        }
    }
}

double scalar(const double* a, const double* b, int N) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += a[i] * b[i];
    }
    return res;
}

void fill(double* x, int N, double value) {
    for (int i = 0; i < N; i++) {
        x[i] = value;
    }
}

void dump(double* x, int N) {
    if (N < 10) {
        for (int i = 0; i < N; i++) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }
}

void FillInitialValues(double* A, double* b, int N) {
    double* u = new double[N];
    fill(A, N*N, 1.0);
    for (int i = 0; i < N; i++) {
        A[i * N + i] = 2.0;
    }
    for (int i = 0; i < N; i++) {
        u[i] = sin(2*M_PI*i / N) * 100;
    }
    mulMV(A, u, N, b);
    dump(u, N);
    delete[] u;
}

void solve(const double* A, const double* b, int N, double* x) {
    std::cout << "STARTED SOLVING " << std::endl;
    double* r = new double[N];
    double* z = new double[N];
    double* Az = new double[N];
    double* buf1 = new double[N];

    fill(x, N, 0.0);
    mulMV(A, x, N, buf1);
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - buf1[i];
        z[i] = r[i];
    }
    double r_square = scalar(r, r, N);
    double b_square = scalar(b, b, N);

    int iter;
    for (iter = 0; iter < 100; iter++) {
        std::cout << "ITER: " << iter << std::endl;

        mulMV(A, z, N, Az);

        double alpha = r_square / scalar(Az, z, N);
        for (int i = 0; i < N; i++) {
            x[i] = x[i] + alpha * z[i];
        }

        for (int i = 0; i < N; i++) {
            buf1[i] = r[i] - alpha * Az[i];
        }

        double r_new_square = scalar(buf1, buf1, N);
        double beta = r_new_square / r_square;

        for (int i = 0; i < N; i++) {
            r[i] = buf1[i];
        }

        r_square = r_new_square;

        for (int i = 0; i < N; i++) {
            z[i] = r[i] + beta * z[i];
        }

        if ((r_square / b_square) < EPS) {
            break;
        }
    }

    std::cout << "DONE" << std::endl;
    delete[] r;
    delete[] z;
    delete[] Az;
    delete[] buf1;
}

int main(int argc, char** argv) {
    const int N = 5;

    // === Initialize MPI ===

    int p_rank;
    int p_count;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);

    const int MAX_PART_SIZE = N / p_count + 1;

    // === Load data and split data across processes ===
    
    double* A_part = new double[MAX_PART_SIZE * N];
    double* b_part = new double[MAX_PART_SIZE];

    double* A;
    double* b;
    double* res;

    if (p_rank == 0) {
        A = new double[N*N];
        b = new double[N];
        res = new double[N];
        FillInitialValues(A, b, N);
    }

    int* b_starts = new int[p_count];
    int* b_sizes = new int[p_count];
    int* A_starts = new int[p_count];
    int* A_sizes = new int[p_count];

    if (p_rank == 0) {
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
    }

    std::cout << "process " << p_rank << " will do " << b_starts[p_rank] << "-" << b_starts[p_rank]+b_sizes[p_rank]  << std::endl;

    MPI_Scatterv(A, A_sizes, A_starts, MPI_DOUBLE, A_part, A_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b, b_sizes, b_starts, MPI_DOUBLE, b_part, b_sizes[p_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // !TODO: change to "by parts"
    solve(A, b, N, res);

    std::cout << "RESULT" << std::endl;
    dump(res, N);

    if (p_rank == 0) {
        delete[] A;
        delete[] b;
        delete[] res;
    }

    delete[] b_starts;
    delete[] b_sizes;
    delete[] A_starts;
    delete[] A_sizes;

    return 0;
}