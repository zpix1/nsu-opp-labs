import os
import numpy as np
import matplotlib.pyplot as plt

data = {'2_T': [72.3259, 36.5395, 25.1762, 12.6711, 15.9135, 13.2452], '2_Sp': [1.0, 1.9793894278794184, 2.872788586045551, 5.7079416940912795, 4.544939830961134, 5.460536647238245], '1_PC': [1, 2, 4, 8, 12, 16, 24], '2_PC': [1, 2, 4, 8, 12, 16], '1_Ep': [100.0, 64.90253673960692, 99.11179210859866, 46.91572946157243, 35.878417164729036, 36.367237519227515, 38.58482968697265], '2_Ep': [100.0, 98.96947139397092, 71.81971465113878, 71.34927117614099, 37.87449859134278, 34.12835404523903], '1_Sp': [1.0, 1.2980507347921382, 3.9644716843439465, 3.7532583569257945, 4.3054100597674845, 5.818758003076403, 9.260359124873435], '1_T': [72.2521, 55.662, 18.2249, 19.2505, 16.7817, 12.4171, 7.8023]}

def gen_plot(x, y, plot, name, v):
    if name == 'T':
        plot.set_title('Время выполнения, с.' + ' (в. {})'.format(v))
        plot.set_ylim([0, 70])
    if name == 'Sp':
        plot.set_title('Ускорение' + ' (в. {})'.format(v))
        plot.set_ylim([0, 10])
    if name == 'Ep':
        plot.set_title('Эффективность, %' + ' (в. {})'.format(v))
        plot.set_ylim([0, 100])
    plot.plot(x, y)
    plot.set_xticks(x)
    return plot

fig, plots = plt.subplots(nrows=2, ncols=3)

print(plots)

i = 0
for v in ['1', '2']:
    for t in ['T', 'Sp', 'Ep']:

        x = data['{}_{}'.format(v, 'PC')]
        y = data['{}_{}'.format(v, t)]
        gen_plot(x, y, plots.flatten()[i], t, v)
        i += 1
fig.tight_layout()
fig.savefig('results_plots.png', dpi=900)