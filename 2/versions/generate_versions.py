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

types = [1, 2]
cores = [1, 2, 4, 8, 12, 16, 24]

def write_version(nodes, cpus, version):
    npc = nodes * cpus
    script = base.format(NODES=nodes, CPUS=min(cpus, 24), NPC=npc, VERSION=version)
    fname = f'run_mul_version_{version}_{npc}_{nodes}-{cpus}.sh'
    with open(fname, 'w') as f:
        f.write(script)
    os.system(f'chmod +x {fname}')
    print(fname)

write_version(1, 1, 1)
write_version(1, 2, 1)
write_version(1, 4, 1)
write_version(1, 8, 1)
write_version(1, 12, 1)
write_version(1, 16, 1)
write_version(1, 20, 1)
write_version(1, 24, 1)
write_version(1, 28, 1)
write_version(1, 32, 1)

write_version(1, 1, 2)
write_version(1, 2, 2)
write_version(1, 4, 2)
write_version(1, 8, 2)
write_version(1, 12, 2)
write_version(1, 16, 2)
write_version(1, 20, 2)
write_version(1, 24, 2)
write_version(1, 28, 2)
write_version(1, 32, 2)