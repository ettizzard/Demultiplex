#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=r2plot     ### Job Name
#SBATCH --output=r2plot.out        ### File in which to store job output
#SBATCH --error=r2plot.err          ### File in which to store job error messages
#SBATCH --time=0-08:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --cpus-per-task=8
#SBATCH --account=bgmp      ### Account used for job submission
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tizzard@uoregon.edu
 
set -e

/usr/bin/time -v ./perbaseqsdistributions.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -l 8 -r 2