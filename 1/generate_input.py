import os
import numpy as np
import matplotlib.pyplot as plt

from settings import *

def random_matrix():
    A_t = np.random.uniform(0, 10, size=(N, N))
    A = (A_t + A_t.T)/2

    x = np.random.uniform(0, 10, size=(N,))

    b = A.dot(x)

    return A, b, x

def plate_matrix():
    P = np.zeros(shape=(Nx, Ny))
    Ox = Nx // 2
    Oy = Ny // 2
    R = 10
    for i in range(Nx):
        for j in range(Ny):
            y = i - Ox
            x = j - Oy

            x /= 10
            y /= -10
            
            if i == 0 or j == 0 or i == Nx - 1 or j == Ny - 1:
                P[i][j] = 1
            elif (R) ** 2 <= (i - Ox)**2 + (j - Oy)**2 <= (R+1) ** 2:
                P[i][j] = -1
            # if (x**2 + y**2 - 1)**3 - x**2*y**3 < 0:
            #     P[i][j] = -1
    # P[Ox][Oy] = 10
    return P

def neighbors(i, j):
    ans = []
    nb = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    for n in nb:
        ni = i + n[0]
        nj = j + n[1]
        if ni < Nx and nj < Ny:
            if ni >= 0 and nj >= 0:
                ans.append((ni, nj))
    return ans

def plate_sol_matrix():
    A = np.zeros(shape=(N, N))
    for i in range(Nx):
        for j in range(Ny):
            idx = i * Ny + j
            nb = neighbors(i, j)
            A[idx][idx] = -4
            for n in nb:
                nidx = n[0] * Ny + n[1]
                A[idx][nidx] = 1
    return A

P = plate_matrix()

A = plate_sol_matrix()
b = P.flatten()

# print(A)
# print(A == A.T)
# print(b)

with open('input.dat', 'w') as f:
    f.write(str(b.size) + '\n')
    A.tofile(f, sep=" ")
    f.write('\n')
    b.tofile(f, sep=" ")

