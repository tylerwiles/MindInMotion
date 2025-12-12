#!/bin/bash

# Required Items
#SBATCH --nodes=1 # Number of Physical Servers
#SBATCH --ntasks=1 # Number of CPUs
#SBATCH --cpus-per-task=24 # Number of cores per task
#SBATCH --mem=256gb # Job memory request
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically 
#SBATCH --time=0-06:00 # Time limit days-hrs:min
#SBATCH --qos=dferris # -b to use Burst

# Optional Items
#SBATCH --job-name=HSIM_Gaitprint # Job name 
#SBATCH --mail-type=END,FAIL # Mail events 
#SBATCH --mail-user=twiles@ufl.edu # Where to send mail 
#SBATCH --output=HSIM_Gaitprint%j.log # Standard output and error log 
pwd; hostname; date # Some user information

# Run the code
module purge
module load R # Load needed modules

# Important for packages
export R_LIBS_USER="/blue/dferris/twiles/Rlibs"
mkdir -p "$R_LIBS_USER"
echo "Using personal R library at: $R_LIBS_USER"

echo "Running H Simulation Test" # Tell me its doing something

Rscript /blue/dferris/twiles/MindInMotion/R/H_Simulation_Gaitprint.R

date # Print ending time
