#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <mpi.h>
#include <mpe.h>

const double EPS = 10e-8;

const int K = 3;

const int Nx = 144*K + 2;
const int Ny = 144*K + 2;
const int Nz = 144*K + 2;

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

inline double rho(int i, int j, int k) {
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


    double denom = 2 * (1 / (Hx*Hx) + 1 / (Hy*Hy) + 1 / (Hz*Hz)) + a;

    double numer =  (phi_i_p1 + phi_i_m1) / (Hx * Hx) +
                    (phi_j_p1 + phi_j_m1) / (Hy * Hy) +
                    (phi_k_p1 + phi_k_m1) / (Hz * Hz) - rho(i, j, k);

    double res = numer / denom;

    // if (res < 0.01) {
    //     DEBUG_PRINT("res is small");
        // DEBUG(res);
        // DEBUG(phi_i_p1);
        // DEBUG(phi_j_p1);
        // DEBUG(phi_k_p1);
        // DEBUG(phi_i_m1);
        // DEBUG(phi_j_m1);
        // DEBUG(phi_k_m1);
    // }

    return res;
}

void print_ans(double* ans, int k0, int k1) {
    if (p_rank == 0) {
        for (int k = 0; k < k1 - k0; k++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                    double xi = x0 + i * Hx;
                    double yi = y0t + j * Hy;
                    double zi = z0 + (k + k0) * Hz;

                    double got = ans[k * Nx * Ny + i + j * Nx];
                    double should = phi(i, j, k+k0);
                    double err =  fabs(got-should);
                    if (err > 10e-5 || Nx < 10) {
                        printf("%9.2lg, %9.2lg, %9.2lg:\t |%9.2lg - %9.2lg| = %lg\n", xi, yi, zi, got, should, fabs(got-should));
                    }
                }
            }
        }
    }
}

// ans should be Nx * Ny * Nz length
void solve(double* ans) {
    const int STEP = (Nz-2) / p_count;
    // z coordinate
    const int BASE = STEP * p_rank;
    const int ZSTRIDE = Nx * Ny;

    double* toplayer = new double[ZSTRIDE]();
    double* bottomlayer = new double[ZSTRIDE]();

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            if (p_rank == p_count - 1)
                toplayer[x + Nx*y] = phi(x, y, Nz - 1);
            else
                toplayer[x + Nx*y] = 0;
            
            if (p_rank == 0)
                bottomlayer[x + Nx*y] = phi(x, y, 0);
            else
                bottomlayer[x + Nx*y] = 0;
        }
    }

    double* baselayer = new double[STEP * ZSTRIDE]();
    double* newbaselayer = new double[STEP * ZSTRIDE]();

    int below_proc = p_rank - 1;
    int above_proc = p_rank + 1 == p_count ? -1 : p_rank + 1;

    // std::cout << "I am " << p_rank << "; Upper me " << above_proc << "; Down me " << below_proc << std::endl;
    
    int stop_iteration = 0;
    for (int i = 0; !stop_iteration; i++) {
        // if (p_rank == 0)
        //     DEBUG(i);

        MPI_Request sent_top_req, sent_bottom_req, got_top_req, got_bottom_req;

    
        // for base
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                newbaselayer[x + y * Nx] = jacobi(  x, y, BASE + 1, 
                                                    bottomlayer, 
                                                    baselayer, 
                                                    STEP > 1 ? baselayer + ZSTRIDE : toplayer);
            }
        }

        // for base + STEP - 1
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                newbaselayer[(STEP - 1) * ZSTRIDE + x + y * Nx] = jacobi(   x, y, BASE + STEP - 0, 
                                                                            STEP > 1 ? baselayer + (STEP - 2) * ZSTRIDE : bottomlayer, 
                                                                            baselayer + (STEP - 1) * ZSTRIDE, 
                                                                            toplayer);
            }
        }

        // async send base and base + STEP - 1
        if (below_proc != -1) {
            MPI_Irecv(bottomlayer, ZSTRIDE, MPI_DOUBLE, below_proc, 1, MPI_COMM_WORLD, &got_bottom_req);
            MPI_Isend(newbaselayer, ZSTRIDE, MPI_DOUBLE, below_proc, 0, MPI_COMM_WORLD, &sent_bottom_req);
        }
        if (above_proc != -1) {
            MPI_Irecv(toplayer, ZSTRIDE, MPI_DOUBLE, above_proc, 0, MPI_COMM_WORLD, &got_top_req);
            MPI_Isend(newbaselayer + (STEP - 1) * ZSTRIDE, ZSTRIDE, MPI_DOUBLE, above_proc, 1, MPI_COMM_WORLD, &sent_top_req);
        }

        // calc center
        for (int dz = 1; dz < STEP - 1; dz++) {
            for (int x = 0; x < Nx; x++) {
                for (int y = 0; y < Ny; y++) {
                    newbaselayer[dz * ZSTRIDE + x + y * Nx] = jacobi(   x, y, BASE + dz + 1, 
                                                                        baselayer + (dz - 1) * ZSTRIDE, 
                                                                        baselayer + (dz + 0) * ZSTRIDE, 
                                                                        baselayer + (dz + 1) * ZSTRIDE);
                }
            }
        }


        int converges = 1;

        for (int i = 0; i < STEP * ZSTRIDE; i++) {
            if (fabs(baselayer[i] - newbaselayer[i]) > EPS) {
                converges = 0;
                break;
            }
        }
        
        if (below_proc != -1) {
            MPI_Wait(&sent_bottom_req, MPI_STATUS_IGNORE);
        }
        if (above_proc != -1) {
            MPI_Wait(&sent_top_req, MPI_STATUS_IGNORE);
        }

        std::swap(newbaselayer, baselayer);

        MPI_Allreduce(&converges, &stop_iteration, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

        if (below_proc != -1) {
            MPI_Wait(&got_bottom_req, MPI_STATUS_IGNORE);
        }
        if (above_proc != -1) {
            MPI_Wait(&got_top_req, MPI_STATUS_IGNORE);
        }
    }

    #ifdef CHECK
    if (p_rank == 0) {
        print_ans(bottomlayer, 0, 1);
        print_ans(baselayer, 1, STEP + 1);
        print_ans(toplayer, STEP + 1, STEP + 2);
    }
    #endif
    
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
    MPE_Init_log();
    
    double begin, end;

    double* ans;

    if (p_rank == 0) {
        begin = MPI_Wtime();
        ans = new double[Nx * Ny * Nz];
    }

    solve(ans);
    
    if (p_rank == 0) {
        end = MPI_Wtime();
        std::cout << p_count << "," << end - begin << std::endl;
        delete[] ans;
    }
    MPE_Finish_log("main1mpe");
    MPI_Finalize();
}