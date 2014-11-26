#!/bin/bash
#PBS -l select=1:ncpus=2
#PBS -l walltime=00:01:00

# to run: qsub 'script_name'.sh
# ensure relevant modules are loaded e.g. mpi, intel-suite ...

# find pglt_farm.py in $PATH
FARM = `which pglt_farm.py`
# run FARM with a directory
mpiexec python $FARM abs_path_to_folder
