#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

const double EPS = 10e-8;

const int Nx = 4;
const int Ny = 4;
const int Nz = 4;

const double Dx = 2;
const double Dy = 2;
const double Dz = 2;

const double x0 = -1.;
const double y0t = -1.;
const double z0 = -1.;

const double Hx = Dx / (Nx - 1);
const double Hy = Dy / (Ny - 1);
const double Hz = Dz / (Nz - 1);

const double a = 10e5;

int p_rank;
int p_count;

#define DEBUG(var) \
            do { std::cout << #var << ": " << var << std::endl; } while (0)

#define DEBUG_PRINT(str) \
            do { std::cout << str << std::endl; } while (0)

inline double phi(int i, int j, int k) {
    double xi = x0 + i * Hx;
    double yi = y0t + j * Hy;
    double zi = z0 + k * Hz;

    return xi*xi + yi*yi + zi*zi;
}

inline double ro(int i, int j, int k) {
    return 6 - a * phi(i, j, k);
}

inline double jacobi(int i, int j, int k, double* bottom, double* cur, double* next) {
    if (((i == 0 || i == Nx - 1) || (j == 0 || j == Ny - 1)) || (k == 0 || k == Nz - 1)) {
        return phi(i, j, k);
    }

    double phi_i_p1 = cur[i+1 + j * Nx];
    double phi_i_m1 = cur[i-1 + j * Nx];

    double phi_j_p1 = cur[i + (j+1) * Nx];
    double phi_j_m1 = cur[i + (j-1) * Nx];

    double phi_k_p1 = next[i + j * Nx];
    double phi_k_m1 = bottom[i + j * Nx];


    double down = 2 * (1 / Hx*Hx + 1 / Hy*Hy + 1 / Hz*Hz) + a;

    double up = (phi_i_p1 + phi_i_m1) / (Hx * Hx) + \
                (phi_j_p1 + phi_j_m1) / (Hy * Hy) + \
                (phi_k_p1 + phi_k_m1) / (Hz * Hz) - ro(i, j, k);
    double res = up / down;

    // if (res < 0.01) {
    //     DEBUG_PRINT("res is small");
        DEBUG(res);
    // }

    return res;
}

void print_ans(double* ans) {
    if (p_rank == 0) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < Nz; k++) {
                    double xi = x0 + i * Hx;
                    double yi = y0t + j * Hy;
                    double zi = z0 + k * Hz;

                    double got = ans[k * Nx * Ny + i + j * Nx];
                    double should = phi(i, j, k);

                    printf("%.2lf, %.2lf, %.2lf:\t |%lf - %lf| = %lf\n", xi, yi, zi, got, should, fabs(got-should));
                }
            }
        }
    }
}

// ans should be Nx * Ny * Nz length
void solve(double* ans) {
    const int STEP = Nz / p_count;
    // z coordinate
    const int BASE = STEP * p_rank;
    const int ZSTRIDE = Nx * Ny;

    double* toplayer = new double[ZSTRIDE];
    double* bottomlayer = new double[ZSTRIDE];

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            toplayer[x + Nx*y] = phi(x, y, BASE + STEP - 1);
            bottomlayer[x + Nx*y] = phi(x, y, BASE);
        }
    }

    double* baselayer = new double[STEP * ZSTRIDE];
    double* newbaselayer = new double[STEP * ZSTRIDE];

    int recv_from = p_rank - 1;
    int send_to = p_rank + 1 == p_count ? -1 : p_rank + 1;

    std::cout << "I am " << p_rank << "; I send to " << send_to << "; I recv from " << recv_from << std::endl;
    
    int flag = 1;
    for (int i = 0; flag; i++) {
        if (p_rank == 0)
            DEBUG(i);

        MPI_Request req_bottom, req_top;

        // for base
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                newbaselayer[x + y * Nx] = jacobi(  x, y, BASE, 
                                                    bottomlayer, 
                                                    baselayer, 
                                                    baselayer + ZSTRIDE);
            }
        }

        // for base + STEP - 1
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                newbaselayer[(STEP - 1) * ZSTRIDE + x + y * Nx] = jacobi(   x, y, BASE + STEP - 1, 
                                                                            baselayer + (STEP - 2) * ZSTRIDE, 
                                                                            baselayer + (STEP - 1) * ZSTRIDE, 
                                                                            toplayer);
            }
        }


        // async send base and base + STEP - 1
        if (send_to != -1) {
            MPI_Isend(newbaselayer, ZSTRIDE, MPI_DOUBLE, send_to, 0, MPI_COMM_WORLD, &req_bottom);
            MPI_Isend(newbaselayer + (STEP - 1) * ZSTRIDE, ZSTRIDE, MPI_DOUBLE, send_to, 1, MPI_COMM_WORLD, &req_top);
        }

        // calc center
        for (int z = BASE + 1; z < BASE - 1; z++) {
            for (int x = 0; x < Nx; x++) {
                for (int y = 0; y < Ny; y++) {
                    newbaselayer[(z - BASE) * ZSTRIDE + x + y * Nx] = jacobi(   x, y, z, 
                                                                                baselayer + (z - BASE - 1) * ZSTRIDE, 
                                                                                baselayer + (z - BASE) * ZSTRIDE, 
                                                                                baselayer + (z - BASE + 1) * ZSTRIDE);
                }
            }
        }

        double max = 0;
        if (p_rank == 0) {
            print_ans(newbaselayer);
            for (int i = 0; i < STEP * ZSTRIDE; i++) {
                if (fabs(baselayer[i] - newbaselayer[i]) > max) {
                    max = fabs(baselayer[i] - newbaselayer[i]);
                }
                baselayer[i] = newbaselayer[i];
            }
            if (max < EPS) {
                flag = 0;
            }
        } else {
            for (int i = 0; i < STEP * ZSTRIDE; i++) {
                baselayer[i] = newbaselayer[i];
            }
        }

        MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // receive next toplayer and next bottomlayer
        if (recv_from != -1) {
            MPI_Irecv(toplayer, ZSTRIDE, MPI_DOUBLE, recv_from, 0, MPI_COMM_WORLD, &req_bottom);
            MPI_Irecv(bottomlayer, ZSTRIDE, MPI_DOUBLE, recv_from, 1, MPI_COMM_WORLD, &req_top);

            MPI_Wait(&req_bottom, MPI_STATUS_IGNORE);
            MPI_Wait(&req_top, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_Gather(newbaselayer, STEP*ZSTRIDE, MPI_DOUBLE, ans, STEP*ZSTRIDE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] toplayer;
    delete[] bottomlayer;
    delete[] baselayer;
}

int main(int argc, char** argv) {
    // INIT MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    
    double begin, end;

    double* ans;

    if (p_rank == 0) {
        begin = MPI_Wtime();
        ans = new double[Nx * Ny * Nz];
    }

    solve(ans);
    
    if (p_rank == 0) {
        end = MPI_Wtime();
        std::cout << "DONE: " << end - begin << std::endl;
        print_ans(ans);
        delete[] ans;
    }

    MPI_Finalize();
}