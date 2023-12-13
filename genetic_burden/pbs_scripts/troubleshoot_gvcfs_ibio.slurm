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

#gatk CombineGVCFs -R /home/durwa004/durwa004/GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fasta --variant M6468.gatk.gvcf.gz --variant M989.gatk.gvcf.gz -O ibio.g.vcf.gz

gatk --java-options "-Xmx4g" GenotypeGVCFs -R /home/durwa004/durwa004/GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fasta -V ibio.g.vcf.gz -O ibio.vcf.gz

bcftools stats ibio.vcf.gz > ibio.stats
bcftools stats -s M6468,M989  ibio.vcf.gz > ind_ibio.stats
