#!/bin/bash
#PBS -l select=1:ncpus=12
#PBS -l walltime=00:01:00

# ensure relevant modules are loaded e.g. mpi, intel-suite ...
# to run: qsub -v folder='abs_path_to_parent_folder' 'script_name'.sh

# find pglt_farm.py in $PATH
FARMSCRIPT=`which pglt_farm.py`
# run $FARMSCRIPT with an absolute path to a parent folder
FOLDER=${folder}
echo 'Running with' $FOLDER
mpiexec python $FARMSCRIPT $FOLDER
