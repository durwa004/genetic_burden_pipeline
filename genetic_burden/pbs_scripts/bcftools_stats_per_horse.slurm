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

#bcftools stats -S horse_IDs.txt /panfs/jay/groups/6/durwa004/shared/PopulationVCF/joint_genotype.goldenPath.vep.vcf.gz > bcftools_stats_per_horse.stats

#bcftools stats -S horse_IDs.txt /panfs/jay/groups/6/durwa004/shared/PopulationVCF/joint_genotype_indels.goldenPath.vep.vcf.gz > bcftools_stats_indels_per_horse.stats

bcftools stats -S horse_IDs.txt /panfs/jay/groups/6/durwa004/shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz > bcftools_stats_combined_per_horse.stats
