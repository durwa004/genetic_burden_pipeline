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

conda activate ensembl-vep
module load bcftools
cd ~/genetic_burden/ 
bcftools stats -S horse_IDs.txt joint_genotype_subset.vep.vcf.gz > joint_genotype_subset.vep.ind.stats

cd /scratch.global/marlo072/CheckHorses/thesis_union/
bcftools stats -S /home/durwa004/durwa004/genetic_burden/horse_IDs.txt /scratch.global/marlo072/CheckHorses/thesis_union/thesis_union.vcf.gz > /scratch.global/marlo072/CheckHorses/thesis_union/thesis_union.ind.stats

