#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>


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
#define DEBUG(var) \
            do { std::cout << #var << ": " << var << std::endl; } while (0)

void dump(double* x, int N) {
    if (N < 20) {
        for (int i = 0; i < N; i++) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }
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
    DEBUG(r_square);
    DEBUG(b_square);
    DEBUG("R");
    dump(r, N);


    int iter;
    for (iter = 0; iter < 100000; iter++) {
        std::cout << "ITER: " << iter << std::endl;

        mulMV(A, z, N, Az);
        dump(Az, N);

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

        for (int i = 0; i < N; i++) {
            z[i] = r[i] + beta * z[i];
        }

        if ((sqrt(r_square) / sqrt(b_square)) < EPS) {
            break;
        }

        r_square = r_new_square;
    }

    std::cout << "DONE" << std::endl;
    delete[] r;
    delete[] z;
    delete[] Az;
    delete[] buf1;
}

int main() {
    const int N = 1000;
    // === 1 ===
    // Matrix A(N);
    // A.init(1.0);
    // for (int i = 0; i < N; i++) {
    //     A.set(i, i, 2.0);
    // }
    // VD b(N, N+1);
    // VD res = solve(A, b);

    // === 2 ===
    double* A = new double[N*N];
    double* u = new double[N];
    double* b = new double[N];
    double* res = new double[N];
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
    solve(A, b, N, res);

    // std::cout << "RESULT" << std::endl;
    // for (int i = 0; i < N; i++) {
    //     std::cout << res[i] << " ";
    // }
    // std::cout << std::endl;

    // mulMV(A, res, N, u);
    // for (int i = 0; i < N; i++) {
    //     std::cout << fabs(u[i] - b[i]) << std::endl;
    // }

    delete[] A;
    delete[] u;
    delete[] b;
    delete[] res;

    return 0;
}