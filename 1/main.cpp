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

int main() {
    const int N = 5;
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
    
    fill(A, N*N, 1.0);
    for (int i = 0; i < N; i++) {
        A[i * N + i] = 2.0;
    }
    for (int i = 0; i < N; i++) {
        u[i] = sin(2*M_PI*i / N) * 100;
    }
    mulMV(A, u, N, b);
    solve(A, b, N, res);
    std::cout << "VECTOR u: " << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "RESULT" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << res[i] << " ";
    }
    std::cout << std::endl;

    delete[] A;
    delete[] u;
    delete[] b;
    delete[] res;

    return 0;
}