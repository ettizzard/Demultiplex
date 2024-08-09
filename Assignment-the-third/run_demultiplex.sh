#!/bin/bash

#SBATCH --account=bgmp   ### change this to your actual account for charging
#SBATCH --partition=bgmp     ### queue to submit to
#SBATCH --job-name=demultiplex    ### job name
#SBATCH --output=demultiplex.out   ### file in which to store job stdout
#SBATCH --error=demultiplex.err    ### file in which to store job stderr
#SBATCH --time=1-0                ### wall-clock time limit, in minutes
#SBATCH --nodes=1               ### number of nodes to use
#SBATCH --cpus-per-task=1       ### number of cores for each task

./demultiplex.py -i /projects/bgmp/shared/2017_sequencing/indexes.txt \
    -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
    -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
    -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
    -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz