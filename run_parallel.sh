#!/bin/bash
#PBS -l select=1:ncpus=2:mem=2000mb
#PBS -l walltime=02:00:00
#PBS -N PREDICTS_PD_ANALYSIS_TEST

# to run: qsub run_parallel.sh

module load mpi
mpiexec python task_farm.py
