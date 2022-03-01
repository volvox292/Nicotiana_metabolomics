#! /bin/bash
#SBATCH -N 1 --exclusive

module load singularity

# The command file with 1 command per line
jobs=commands.txt

# This script must be launch via "SBATCH --array=1-48" 
hpc_multilauncher $jobs