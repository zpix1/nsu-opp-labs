#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <omp.h>
const double EPS = 1e-5;

void dump(double* x, int N, std::string fname) {
    std::ofstream f(fname);
    for (int i = 0; i < N; i++) {
        f << x[i] << " ";
    }
    f << std::endl;
    f.close();
}

// void mulMV(const double* A, const double* x, int N, double* res)  {
//     for (int i = 0; i < N; i++) {
//         res[i] = 0;
//         for (int j = 0; j < N; j++) {
//             res[i] += A[i * N + j] * x[j];
//         }
//     }
// }

// double scalar(const double* a, const double* b, int N) {
//     double res = 0;
//     for (int i = 0; i < N; i++) {
//         res += a[i] * b[i];
//     }
//     return res;
// }

void fill(double* x, int N, double value) {
    for (int i = 0; i < N; i++) {
        x[i] = value;
    }
}

void solve(const double* A, const double* b, int N, double* x) {
    double* r = new double[N];
    double* z = new double[N];
    double* Az = new double[N];
    double* buf1 = new double[N];

    int i, j;

    for (i = 0; i < N; i++) {x[i] = 0.0;}
    // mulMV(A, x, N, buf1);
    for (i = 0; i < N; i++) {
        buf1[i] = 0;
        for (j = 0; j < N; j++) {
            buf1[i] += A[i * N + j] * x[j];
        }
    }

    for (i = 0; i < N; i++) {
        r[i] = b[i] - buf1[i];
        z[i] = r[i];
    }

    // double r_square = scalar(r, r, N);
    double r_square = 0.0;
    for (i = 0; i < N; i++) { r_square += r[i] * r[i]; }
    // double b_square = scalar(b, b, N);
    double b_square = 0.0;
    for (i = 0; i < N; i++) { b_square += b[i] * b[i]; }

    // int iter;
    double sc_Az_z = 0, r_new_square = 0;
    #pragma omp parallel shared(sc_Az_z, r_new_square)
    for (int iter = 0; iter < 100000; iter++) {
        #pragma omp single
        std::cout << "ITER: " << iter << std::endl;
            // mulMV(A, z, N, Az);
        for (i = 0; i < N; i++) {
            Az[i] = 0;
            for (j = 0; j < N; j++) {
                Az[i] += A[i * N + j] * z[j];
            }
        }

        sc_Az_z = 0;
        #pragma omp barrier

        double sc_Az_z_private = 0;
        #pragma omp for 
        for (i = 0; i < N; i++) { sc_Az_z_private += Az[i] * z[i]; }

        #pragma omp atomic
        sc_Az_z += sc_Az_z_private;
        #pragma omp barrier
        
        double alpha = r_square / sc_Az_z;
        
        #pragma omp for
        for (i = 0; i < N; i++) {
            x[i] = x[i] + alpha * z[i];
            buf1[i] = r[i] - alpha * Az[i];
        }
        r_new_square = 0;
        #pragma omp barrier

        // double r_new_square = scalar(buf1, buf1, N);
        double r_new_square_private = 0;
        #pragma omp for
        for (i = 0; i < N; i++) { r_new_square_private += buf1[i] * buf1[i]; }

        #pragma omp atomic
        r_new_square += r_new_square_private;
        #pragma omp barrier

        double beta = r_new_square / r_square;
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < N; i++) {
            r[i] = buf1[i];
        }
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < N; i++) {
            z[i] = r[i] + beta * z[i];
        }
        if ((sqrt(r_square) / sqrt(b_square)) < EPS) {
            break;
        }
        #pragma omp barrier

        r_square = r_new_square;
        #pragma omp barrier
    }

    delete[] r;
    delete[] z;
    delete[] Az;
    delete[] buf1;
}

int main() {
    
    int N;
    
    std::ifstream f("input.dat");
    f >> N;
    double* A = new double[N*N];
    double* b = new double[N];
    double* res = new double[N];
    for (int i = 0; i < N*N; i++) {
        f >> A[i];
    }
    for (int i = 0; i < N; i++) {
        f >> b[i];
    }

    double start = omp_get_wtime();
    solve(A, b, N, res);
    double end = omp_get_wtime();
    std::cout << "DONE: " << (end - start) << std::endl;

    dump(res, N, "output.dat");

    delete[] A;
    delete[] b;
    delete[] res;
    return 0;
}