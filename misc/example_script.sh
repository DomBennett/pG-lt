#!/bin/bash
#PBS -l select=1:ncpus=12
#PBS -l walltime=00:01:00

# ensure relevant modules are loaded e.g. mpi, intel-suite ...
# to run: qsub 'script_name'.sh 'abs_path_to_parent_folder'

# find pglt_farm.py in $PATH
FARM=`which pglt_farm.py`
# run FARM with an absolute path to a parent folder
echo 'Running with ' $1
mpiexec python $FARM $1
