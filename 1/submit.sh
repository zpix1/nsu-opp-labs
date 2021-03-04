#!/bin/bash

#PBS -l select=2:ncpus=8:mpiprocs=16:mem=2000m,place=free
#PBS -l walltime=00:01:00
#PBS -q S3077545
# cp runner.py settings.py generate_input.py main*.cpp run.sh $PBS_O_WORKDIR

cd $PBS_O_WORKDIR
echo "I run on node: `uname -n`"
echo "My working directory is: $PBS_O_WORKDIR"
# echo "Assigned to me nodes are:"
# cat $PBS_NODEFILE

echo "Setting conda"
source /opt/shared/anaconda/anaconda3-2020/bin/activate
echo "Generating input"
python3 generate_input.py
echo "Running"
python3 runner.py