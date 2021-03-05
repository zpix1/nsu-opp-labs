#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <math.h>
#include <omp.h>
const double EPS = 1e-10;

void dump(double* x, int N, std::string fname) {
    std::ofstream f(fname);
    for (int i = 0; i < N; i++) {
        f << x[i] << " ";
    }
    f << std::endl;
    f.close();
}

void mulMV(const double* A, const double* x, int N, double* res)  {
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < N; i++) {
        res[i] = 0;
        for (int j = 0; j < N; j++) {
            res[i] += A[i * N + j] * x[j];
        }
    }
}

double scalar(const double* a, const double* b, int N) {
    
    double res = 0;
    #pragma omp parallel for schedule(runtime) reduction(+:res)
    for (int i = 0; i < N; i++) {
        res += a[i] * b[i];
    }
    return res;
}

void fill(double* x, int N, double value) {
    #pragma omp parallel for schedule(runtime)
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

    fill(x, N, 0.0);
    mulMV(A, x, N, buf1);
    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - buf1[i];
        z[i] = r[i];
    }
    double r_square = scalar(r, r, N);
    double b_square = scalar(b, b, N);

    int iter;
    for (iter = 0; iter < 100000; iter++) {
        std::cout << "ITER: " << iter << std::endl;

        mulMV(A, z, N, Az);

        double alpha = r_square / scalar(Az, z, N);

        #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < N; i++) {
            x[i] = x[i] + alpha * z[i];
        }
        #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < N; i++) {
            buf1[i] = r[i] - alpha * Az[i];
        }

        double r_new_square = scalar(buf1, buf1, N);
        double beta = r_new_square / r_square;
        #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < N; i++) {
            r[i] = buf1[i];
        }
        #pragma omp parallel for schedule(runtime)
        for (int i = 0; i < N; i++) {
            z[i] = r[i] + beta * z[i];
        }

        if ((sqrt(r_square) / sqrt(b_square)) < EPS) {
            break;
        }

        r_square = r_new_square;
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