#!/bin/bash
#SBATCH --nodes=1 #request one node
#SBATCH --cpus-per-task=4
#SBATCH --time=2-0:00:00 #ask that the job be allowed to run for 2 days.
#SBATCH --error=job.%J.err # tell it to store the output console text to a file
#SBATCH --output=job.%J.out #tell it to store the error messages to a file
module load R
srun R --vanilla CMD BATCH make2.R

