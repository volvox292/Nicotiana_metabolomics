#!/bin/bash
#SBATCH --partition=fast
#SBATCH --cpus-per-task=15
#SBATCH --mem=260G
#SBATCH --job-name="Sirius"
#SBATCH --output=sir%j.out


#COMMAND LINE TO EXECUTE
srun --exclusive -N 1 -n 1 sirius --maxmz=850 --cores $SLURM_CPUS_PER_TASK -i /home/delser/alltissues15072021-sirius-py.mgf -o /home/delser/sirius_altissues15072021 formula  -c 50 -p orbitrap zodiac fingerid -d BIO canopus