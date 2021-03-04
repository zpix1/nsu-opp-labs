import os
import numpy as np
import matplotlib.pyplot as plt

from settings import *

sol_x = np.fromfile('output.dat', sep=" ")

result_plate = sol_x.reshape((Nx, Ny))

print(result_plate)

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(result_plate, cmap='viridis')
cb = fig.colorbar(cax, ticks=[-5, 0, 5], label='Температура')
cb.ax.set_yticklabels(['Холодно', 'Ноль', 'Горячо']) 
plt.axis('off')
plt.title('Распределение тепла в пластинке')
fig.savefig('plate.png', dpi=900)