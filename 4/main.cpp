#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

const double EPS = 10e-8;

const int Nx = 192;
const int Ny = 192;
const int Nz = 192;

const double Dx = 2;
const double Dy = 2;
const double Dz = 2;

const double a = 10e5;

int p_rank;
int p_count;

#define DEBUG(var) \
            do { std::cout << #var << ": " << var << std::endl; } while (0)

void solve() {
    const int STEP = Nz / p_count;
    const int BASE = STEP * p_rank;
    
    for (int i = 0;; i++) {
        // for base
        for (;;)

        // for base + STEP - 1
        for (;;)

        // async send top and bottop
        // ...

        // calc center
        for (int z = BASE + 1; z < BASE - 1; z++) {
            // ..
        }

        // receive next top and next bottom
        // ...
    }
}

int main(int argc, char** argv) {
    // INIT MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    
    solve();

    MPI_Finalize();
}