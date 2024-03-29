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

const double EPS = 1e-10;

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
void dump(double* x, int N, std::string fname) {
    std::ofstream f(fname);
    for (int i = 0; i < N; i++) {
        f << x[i] << " ";
    }
    f << std::endl;
    f.close();
}

void FillInitialValues(double* A, double* b, int N) {
    srand(5);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            A[i * N + j] = rand() % 100;
            A[j * N + i] = A[i * N + j];
        }
    }
    for (int i = 0; i < N; i++) {
        A[i * N + i] += 500.0;
    }
    for (int i = 0; i < N; i++) {
        b[i] = rand() % 100;
    }
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         std::cout << A[i * N + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
}

int main(int argc, char** argv) {
    double start, end;
    int N;

    // === Initialize MPI ===

    MPI_Init(&argc, &argv);
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

    double b_square_part = 0;
    for (int i = 0; i < b_sizes[p_rank]; i++) {
        part_buf0[i] = 0;
        for (int j = 0; j < N; j++) {
            part_buf0[i] += A_part[i * N + j] * x_part[j];
        }
        part_buf0[i] = b_part[i] - part_buf0[i];
        part_buf1[i] = part_buf0[i];
        b_square_part += b_part[i] * b_part[i];
    }

    MPI_Allreduce(&b_square_part, &b_square, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allgatherv(part_buf0, b_sizes[p_rank], MPI_DOUBLE, r, b_sizes, b_starts, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(part_buf1, b_sizes[p_rank], MPI_DOUBLE, z, b_sizes, b_starts, MPI_DOUBLE, MPI_COMM_WORLD);
    
    r_square = scalar(r, r, N);

    // 2 - main cycle
    // int i = 0;
    int converges = false;
    while (!converges) {
        // if (p_rank == 0) {
        //     std::cout << "ITER " << i << std::endl;
        //     i++;
        // }
        // get Az
        for (int i = 0; i < b_sizes[p_rank]; i++) {
            part_buf0[i] = 0;
            for (int j = 0; j < N; j++) {
                part_buf0[i] += A_part[i * N + j] * z[j];
            }
        }
        MPI_Allgatherv(part_buf0, b_sizes[p_rank], MPI_DOUBLE, Az, b_sizes, b_starts, MPI_DOUBLE, MPI_COMM_WORLD);

        double alpha = r_square / scalar(Az, z, N);
        
        for (int i = 0; i < b_sizes[p_rank]; i++) {
            x_part[i] = x_part[i] + alpha * z[b_starts[p_rank] + i];
        }

        // full_buf - z_n+1
        for (int i = 0; i < N; i++) {
            full_buf[i] = r[i] - alpha * Az[i];
        }
        double r_new_square = scalar(full_buf, full_buf, N);
        double beta = r_new_square / r_square;
        for (int i = 0; i < N; i++) {
            r[i] = full_buf[i];
        }
        
        for (int i = 0; i < N; i++) {
            z[i] = r[i] + beta * z[i];
        }
        if ((sqrt(r_square) / sqrt(b_square)) < EPS) {
            converges = true;
        }
        r_square = r_new_square;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Gatherv(x_part, b_sizes[p_rank], MPI_DOUBLE, x, b_sizes, b_starts, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (p_rank == 0) {
        end = MPI_Wtime();
        std::cout << "DONE: " << p_count << ": " << end - start << std::endl;
        dump(x, N, "output.dat");

        // mulMV(A, x, N, z);
        // for (int i = 0; i < N; i++) {
        //     DEBUG(b[i]);
        //     DEBUG(z[i]);
        //     DEBUG(x[i]);
        //     std::cout << fabs(z[i] - b[i]) << std::endl;
        // }
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