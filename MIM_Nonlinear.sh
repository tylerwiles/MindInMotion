#!/bin/bash

# Required Items
#SBATCH --nodes=1 # Number of Physical Servers
#SBATCH --ntasks=1 # Number of CPUs
#SBATCH --cpus-per-task=20 # Number of cores per task
#SBATCH --mem=32gb # Job memory request
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically 
#SBATCH --time=0-00:30 # Time limit days-hrs:min
#SBATCH --qos=dferris # -b to use Burst

# Optional Items
#SBATCH --job-name=MIM_Nonlinear # Job name 
#SBATCH --mail-type=END,FAIL # Mail events 
#SBATCH --mail-user=twiles@ufl.edu # Where to send mail 
#SBATCH --output=MIM_Nonlinear%j.log # Standard output and error log 
pwd; hostname; date # Some user information

# Run the code
module purge
module load matlab/2024b # Load needed modules

echo "Matlab started up" # Tell me its doing something

srun matlab -nodisplay -nosplash -batch "run('/blue/dferris/twiles/MindInMotion/MATLAB/SI_Nonlinear.m')"

date # Print ending time