#!/bin/bash
#SBATCH --job-name=DimerTSweep                  # replace name
#SBATCH --mail-user=niggeni@zedat.fu-berlin.de  # replace email address
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=20         # memory per cpu-core, more means less gc time
#SBATCH --qos=standard
#SBATCH --time=01:30:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=none
#SBATCH --output=/scratch/niggeni/Dimer.out   # File to which standard Out- will be written
#SBATCH --error=/scratch/niggeni/Dimer.err

module purge
module load Julia
cd /scratch/niggeni/Dimer
export JULIA_NUM_THREADS=32
julia --optimize=3 --math-mode=fast /home/niggeni/JuliaPMFRG/main_Cluster.jl $SLURM_NTASKS
