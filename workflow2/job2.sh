#!/bin/bash
#SBATCH --nodes=1 #request one node
#SBATCH --cpus-per-task=16
#SBATCH --time=200:00:00 #ask that the job be allowed to run for 150 hours.
#SBATCH --error=job.%J.err # tell it to store the output console text to a file
#SBATCH --output=job.%J.out #tell it to store the error messages to a file
module load R
srun R --vanilla CMD BATCH make2.R
make
