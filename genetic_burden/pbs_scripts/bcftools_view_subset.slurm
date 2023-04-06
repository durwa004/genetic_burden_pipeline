#!/bin/bash -l
#SBATCH --job-name=decompose_SNPs.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=12:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=durwa004@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e bcftools_view.%j.subset.e 
#SBATCH -o bcftools_view.%j.sibset.o 


cd ~/genetic_burden/ 
conda activate bcftools
bcftools view -S horse_IDs_ss.txt -c 1 /home/durwa004/shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz > joint_genotype.goldenPath.vep.subset.vcf

conda deactivate bcftools
conda activate GATK
bgzip joint_genotype.goldenPath.vep.subset.vcf
tabix joint_genotype.goldenPath.vep.subset.vcf.gz
bcftools stats joint_genotype.goldenPath.vep.subset.vcf.gz > joint_genotype.goldenPath.vep.subset.stats
bcftools stats -S horse_IDs.txt joint_genotype.goldenPath.vep.subset.vcf.gz > joint_genotype.goldenPath.vep.subset.ind.stats
