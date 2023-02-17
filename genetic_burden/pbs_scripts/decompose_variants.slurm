#!/bin/bash -l
#SBATCH --job-name=decompose_SNPs.slurm
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

cd ~/genetic_burden/ 

gunzip /panfs/jay/groups/6/durwa004/shared/PopulationVCF/joint_genotype.goldenPath.vep.vcf.gz
/home/mccuem/durwa004/bin/vt/vt decompose -s joint_genotype.goldenPath.vep.vcf -o joint_genotype.goldenPath.vep.decomposed.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip joint_genotype.goldenPath.vep.decomposed.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix joint_genotype.goldenPath.vep.decomposed.vcf.gz && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip joint_genotype.goldenPath.vep.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix joint_genotype.goldenPath.vep.vcf.gz
