from subprocess import check_output
import re
import sys
import os

base = '''
#!/bin/bash

# NPC=mpiprocs=select*ncpus
# 1 <= select <= 2

#PBS -l select={NODES}:ncpus={CPUS}:mpiprocs={NPC}:mem=7000m,place=free
#PBS -l walltime=01:01:00

source /opt/intel/mkl/10.0.3.020/tools/environment/mklvarsem64t.sh

cd $PBS_O_WORKDIR

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin

$MPI_PATH/mpicxx main.cpp -o emain{N1}{N2}{N3}{dim1}{dim2}{NODES}{CPUS}{NPC} -std=c++11 -O2 -xhost=native -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
$MPI_PATH/mpirun -machinefile $PBS_NODEFILE -np {NPC} ./emain{N1}{N2}{N3}{dim1}{dim2}{NODES}{CPUS}{NPC} {N1} {N2} {N3} {dim1} {dim2}
'''

def write_version(nodes, cpus, dims, N1, N2, N3):
    npc = nodes * cpus
    script = base.format(NODES=nodes, CPUS=min(cpus, 24), NPC=npc, dim1=dims[0], dim2=dims[1], N1=N1, N2=N2, N3=N3)
    fname = f'rmv_npc{npc}_nodes{nodes}_cpus{cpus}_dims{dims[0]}-{dims[1]}_ns{N1}-{N2}-{N3}.sh'
    with open(fname, 'w') as f:
        f.write(script)
    os.system(f'chmod +x {fname}')
    print(fname)

# 2000 x 2000 - 2000 x 2000
n1 = 6000
n2 = 6000
n3 = 6000

# write_version(1, 1, [1, 1], n1, n2, n3)
# write_version(1, 2, [1, 2], n1, n2, n3)
# write_version(1, 4, [2, 2], n1, n2, n3)
# write_version(1, 8, [2, 4], n1, n2, n3)
# write_version(2, 6, [2, 6], n1, n2, n3)
# write_version(2, 8, [4, 4], n1, n2, n3)
# write_version(2, 12, [4, 6], n1, n2, n3)

# write_version(2, 12, [1, 24], n1, n2, n3)
# write_version(2, 12, [2, 12], n1, n2, n3)
# write_version(2, 12, [3, 8], n1, n2, n3)
# write_version(2, 12, [4, 6], n1, n2, n3)
# write_version(2, 12, [6, 4], n1, n2, n3)
# write_version(2, 12, [8, 3], n1, n2, n3)
# write_version(2, 12, [12, 2], n1, n2, n3)
# write_version(2, 12, [24, 1], n1, n2, n3)

for i in range(1,6):
    write_version(2, 8, [4, 4], i*10*16*5*2, i*160*5*2, i*64*5*2)