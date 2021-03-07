from subprocess import check_output
import re
import sys
import os

base = '''
#!/bin/bash

#PBS -l select={NODES}:ncpus={CPUS}:ompthreads={NPC}:mem=4000m,place=free
#PBS -l walltime=01:01:00
cd $PBS_O_WORKDIR
echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"

g++ -Wall -std=c++11 -fopenmp main{VERSION}.cpp -o main{VERSION}{NODES}{CPUS}{NPC}
OMP_NUM_THREADS={NPC} ./main{VERSION}{NODES}{CPUS}{NPC}
'''

oldbase = '''
#!/bin/bash

# NPC=mpiprocs=select*ncpus
# 1 <= select <= 2

#PBS -l select={NODES}:ncpus={CPUS}:mpiprocs={NPC}:mem=4000m,place=free
#PBS -l walltime=01:01:00

cd $PBS_O_WORKDIR
echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"

MPI_PATH=/opt/intel/impi/5.0.1.035/intel64/bin

$MPI_PATH/mpicxx main{VERSION}.cpp -o main{VERSION}{NODES}{CPUS}{NPC} -std=c++11
$MPI_PATH/mpirun -machinefile $PBS_NODEFILE -np {NPC} ./main{VERSION}{NODES}{CPUS}{NPC}
'''

def write_version(nodes, cpus, version):
    npc = nodes * cpus
    script = base.format(NODES=nodes, CPUS=min(cpus, 24), NPC=npc, VERSION=version)
    fname = f'rmv_{version}_{npc}_{nodes}-{cpus}.sh'
    with open(fname, 'w') as f:
        f.write(script)
    os.system(f'chmod +x {fname}')
    print(fname)

write_version(1, 1, 1)
write_version(1, 2, 1)
write_version(1, 4, 1)
write_version(1, 8, 1)
write_version(2, 6, 1)
write_version(2, 8, 1)
write_version(2, 10, 1)
write_version(2, 12, 1)
# write_version(4, 7, 1)
# write_version(4, 8, 1)

write_version(1, 1, 2)
write_version(1, 2, 2)
write_version(1, 4, 2)
write_version(1, 8, 2)
write_version(2, 6, 2)
write_version(2, 8, 2)
write_version(2, 10, 2)
write_version(2, 12, 2)
# write_version(4, 7, 2)
# write_version(4, 8, 2)