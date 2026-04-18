#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=4:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --job-name=subsamp-test
#SBATCH --mail-user=dcking@colostate.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
#SBATCH --output=%x.%j.log # gives slurm.ID.log

hash nextflow 2>/dev/null || module load nextflow

nextflow run main.nf  -with-trace trace.txt -profile testing  -stub-run -resume 
