#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <omp.h>
const double EPS = 1e-8;

void dump(double* x, int N, std::string fname) {
    std::ofstream f(fname);
    for (int i = 0; i < N; i++) {
        f << x[i] << " ";
    }
    f << std::endl;
    f.close();
}

void fill(double* x, int N, double value) {
    for (int i = 0; i < N; i++) {
        x[i] = value;
    }
}

#define DEBUG(var) \
            do { std::cout << #var << ": " << var << std::endl; } while (0)

void solve(const double* A, const double* b, int N, double* x) {
    double* r = new double[N];
    double* z = new double[N];
    double* Az = new double[N];
    double* buf1 = new double[N];

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {x[i] = 0.0;}
    // mulMV(A, x, N, buf1);
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        buf1[i] = 0;
        for (int j = 0; j < N; j++) {
            buf1[i] += A[i * N + j] * x[j];
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - buf1[i];
        z[i] = r[i];
    }

    // double r_square = scalar(r, r, N);
    double r_square = 0.0;
    #pragma omp parallel for reduction(+:r_square)
    for (int i = 0; i < N; i++) { r_square += r[i] * r[i]; }
    // double b_square = scalar(b, b, N);
    double b_square = 0.0;
    #pragma omp parallel for reduction(+:b_square)
    for (int i = 0; i < N; i++) { b_square += b[i] * b[i]; }

    double sc_Az_z = 0, r_new_square = 0;
    bool converges = false;
    int iter = 0;
    #pragma omp parallel
    while (!converges) {
        #pragma omp single
        {
            sc_Az_z = 0;
            r_new_square = 0;
            std::cout << "ITER: " << iter << std::endl;
            iter++;
        }

        #pragma omp for
        for (int i = 0; i < N; i++) {
            Az[i] = 0;
            for (int j = 0; j < N; j++) {
                Az[i] += A[i * N + j] * z[j];
            }
        }

        #pragma omp for reduction(+:sc_Az_z)
        for (int i = 0; i < N; i++) {
            sc_Az_z += Az[i] * z[i]; 
        }
        
        double alpha = r_square / sc_Az_z;
        
        #pragma omp for
        for (int i = 0; i < N; i++) {
            x[i] = x[i] + alpha * z[i];
            buf1[i] = r[i] - alpha * Az[i];
        }

        #pragma omp for reduction(+:r_new_square)
        for (int i = 0; i < N; i++) { r_new_square += buf1[i] * buf1[i]; }

        double beta = r_new_square / r_square;

        #pragma omp for
        for (int i = 0; i < N; i++) {
            r[i] = buf1[i];
            z[i] = r[i] + beta * z[i];
        }
        
        #pragma omp single
        if ((sqrt(r_square) / sqrt(b_square)) < EPS) {
            converges = true;
            DEBUG(converges);
        }

        #pragma omp single
        {
            r_square = r_new_square;
        }
    }
    DEBUG("end");

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