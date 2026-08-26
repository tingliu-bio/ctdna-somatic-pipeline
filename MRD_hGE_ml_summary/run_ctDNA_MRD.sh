#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J run_ctDNA

work_path=$SLURM_SUBMIT_DIR
cd $work_path

module load R/4.1.0

Rscript ctDNAtool_MRD_cluster_v3.R 
