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

conda activate GATK
module load bcftools
cd /scratch.global/marlo072/CheckHorses/ 

gatk CombineGVCFs -R /panfs/jay/groups/6/durwa004/shared/VEPDocs/goldenPath.Ec_build-3.0_wMSY.fa --variant M6468.goldenPath.g.vcf.gz --variant M989.goldenPath.g.vcf.gz -O goldenPath.g.vcf.gz

gatk --java-options "-Xmx4g" GenotypeGVCFs -R /panfs/jay/groups/6/durwa004/shared/VEPDocs/goldenPath.Ec_build-3.0_wMSY.fa -V goldenPath.g.vcf.gz -O goldenPath.vcf.gz

bcftools stats goldenPath.vcf.gz > goldenPath.vcf.gz.stats
bcftools stats -s M6468, M989 goldenPath.vcf.gz > ind_goldenPath.vcf.gz.stats
