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

#module load java/openjdk-17.0.2
#cd /home/durwa004/durwa004/packages/snpEff/

#java -jar snpEff.jar EquCab3_y ../../../shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz > ../../genetic_burden/joint_genotype_combined.goldenPath.snpeff.vcf
#java -jar snpEff.jar EquCab3_y /scratch.global/marlo072/CheckHorses/thesis_intersect/thesis_intersect.vcf.gz > /scratch.global/marlo072/CheckHorses/thesis_intersect/thesis_intersect.snpeff.vcf
#java -jar snpEff.jar EquCab3_y /scratch.global/marlo072/CheckHorses/thesis_union/thesis_union.vcf.gz > /scratch.global/marlo072/CheckHorses/thesis_union/thesis_union.snpeff.vcf 
#java -jar snpEff.jar EquCab3_y /home/durwa004/durwa004/genetic_burden/joint_genotype.goldenPath.vep.subset.vcf.gz  > /home/durwa004/durwa004/genetic_burden/joint_genotype.goldenPath.vep.subset.snpeff.vcf

conda activate ensembl-vep
cd /home/durwa004/durwa004/genetic_burden/
#bgzip joint_genotype.goldenPath.vep.subset.snpeff.vcf
tabix joint_genotype.goldenPath.vep.subset.snpeff.vcf.gz
