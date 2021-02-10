#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>

using namespace std;

const double EPS = 1e-5;
using VD = vector<double>;

struct Matrix {
    int N;
    double* data;
    Matrix(int _N) {
        N = _N;
        data = (double*)malloc(N * N * sizeof(double));
    }
    ~Matrix() {
        free(data);
    }
    void init(double v) {
        for (int i = 0; i < N*N; i++) {
            data[i] = v;
        }
    }
    double get(int i, int j) const {
        return data[i * N + j];
    }
    void set(int i, int j, double v) {
        data[i * N + j] = v;
    }
    VD operator*(const VD& v) const {
        assert(N == v.size());
        VD res(N, 0.0);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                res[i] += data[i * N + j] * v[j];
            }
        }
        return res;
    }
};

double operator*(const VD& a, const VD& b) {
    assert(a.size() == b.size());
    double res = 0;
    for (int i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}

VD operator*(double k, const VD& b) {
    VD res(b.size());
    for (int i = 0; i < b.size(); i++) {
        res[i] = k * b[i];
    }
    return res;
}

VD operator+(double k, const VD& b) {
    VD res(b.size());
    for (int i = 0; i < b.size(); i++) {
        res[i] = k + b[i];
    }
    return res;
}

VD operator-(const VD& a, const VD& b) {
    assert(a.size() == b.size());
    VD res(a.size());
    for (int i = 0; i < a.size(); i++) {
        res[i] = a[i] - b[i];
    }
    return res;
}

VD operator+(const VD& a, const VD& b) {
    assert(a.size() == b.size());
    VD res(a.size());
    for (int i = 0; i < a.size(); i++) {
        res[i] = a[i] + b[i];
    }
    return res;
}

VD solve(const Matrix& A, const VD& b) {
    cout << "STARTED SOLVING " << endl;
    VD x(A.N, 0.0);
    VD r = b - A*x;
    VD z = r;
    double r_square = r * r;
    double b_square = b * b;
    int iter;
    for (iter = 0; iter < 100; iter++) {
        cout << "ITER: " << iter << endl;
        VD Az = A*z;
        double alpha = r_square / (Az * z);
        x = x + alpha * z;

        VD new_r = r - alpha * Az;
        double r_new_square = new_r * new_r;
        double beta = r_new_square / r_square;
        r = new_r;
        r_square = r_new_square;
        z = r + beta * z;

        if ((r_square / b_square) < EPS)
            break;
    }
    cout << "DONE" << endl;
    return x;
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
    Matrix A(N);
    A.init(1.0);
    for (int i = 0; i < N; i++) {
        A.set(i, i, 2.0);
    }
    VD u(N);
    for (int i = 0; i < N; i++) {
        u[i] = sin(2*M_PI*i / N) * 100;
    }
    VD b = A*u;
    VD res = solve(A, b);
    cout << "VECTOR u: " << endl;
    for (int i = 0; i < N; i++) {
        cout << u[i] << " ";
    }
    cout << endl;

    cout << "RESULT" << endl;
    for (int i = 0; i < N; i++) {
        cout << res[i] << " ";
    }
    cout << endl;

    return 0;
}