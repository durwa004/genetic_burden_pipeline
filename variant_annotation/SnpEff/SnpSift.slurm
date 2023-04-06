#!/bin/bash -l
#SBATCH --job-name=SnpEff.slurm
#SBATCH --nodes=1 #specify one node
#SBATCH --ntasks=1
#SBATCH --time=24:00:00 
#SBATCH --ntasks=2 
#SBATCH --mem=10g 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=durwa004@umn.edu 
#SBATCH -A durwa004 
#SBATCH -e SnpEff.%j.e 
#SBATCH -o SnpEff.%j.o 

module load java/openjdk-13.0.2
cd /home/durwa004/durwa004/packages/snpEff/


#java -Xmx4g -jar SnpSift.jar filter "((ANN[0].IMPACT has 'HIGH') | (ANN[*].IMPACT = 'MODERATE') | (ANN[*].IMPACT = 'LOW'))" ../../genetic_burden/joint_genotype_combined.goldenPath.snpeff.vcf.gz > ../../genetic_burden/joint_genotype_combined.goldenPath.snpeff.hml.vcf
java -Xmx4g -jar SnpSift.jar filter "((ANN[0].IMPACT has 'HIGH') | (ANN[*].IMPACT = 'MODERATE') | (ANN[*].IMPACT = 'LOW'))" ../../genetic_burden/joint_genotype.goldenPath.vep.subset.snpeff.vcf.gz > ../../genetic_burden/joint_genotype.goldenPath.vep.subset.snpeff.hml.vcf


conda activate ensembl-vep
cd /home/durwa004/durwa004/genetic_burden/
bgzip joint_genotype.goldenPath.vep.subset.snpeff.hml.vcf
tabix joint_genotype.goldenPath.vep.subset.snpeff.hml.vcf.gz

conda activate bcftools
bcftools stats joint_genotype.goldenPath.vep.subset.snpeff.hml.vcf.gz > joint_genotype.goldenPath.vep.subset.snpeff.hml.stats
