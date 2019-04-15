#!/bin/bash 
#SBATCH --job-name="threaded"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=100:100:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40000


# Go to the directoy from which our job was launched
# cd $SLURM_SUBMIT_DIR


#run using 24 cores
export OMP_NUM_THREADS=24


# run an application
~/micro_fwi/gitcode/fault_sources/modeling_spatiotemporal.exe
