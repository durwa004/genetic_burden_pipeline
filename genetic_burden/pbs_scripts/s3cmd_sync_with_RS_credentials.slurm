#!/bin/bash -l
#SBATCH --job-name=bcftools_stats_per_horse.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=durwa004@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e bcftools_view.%j.decompose_SNPs.e 
#SBATCH -o bcftools_view.%j.decompose_SNPs.o 

cd /scratch.global/marlo072/

s3cmd  --config /home/mccuem/shared/.local/.s3cfg sync s3://mccue-lab/ibiodatatransfer2019/GCF_002863925.1_EquCab3.0_genomic/ CheckHorses/

