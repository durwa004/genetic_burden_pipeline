# genetic_burden_pipeline
Scripts and tools to estimate the genetic burden as part of my first aim of my thesis.

# bcftools_stats
- Run bcftools for all files in a directory
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/Generate_bcftools_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -e _intersect.vcf.gz
```
- Extract bcftools info
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/Extract_bcftools_stats.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -p bcftools
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_stats.genotyped.vcf.gz.pbs 
```
- Move stats info back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/genotyped_files/joint_bcftools/bcftools_number_of_variants.txt /Users/durwa004/Desktop
```


**Input = Intersect of bcftools/gatk haplotype caller (group calling) on all NC_ chromosomes (exclude unplaced for now).**

NB - will need to re-estimate for the individual calls as well.

**All python scripts that run from command line need source activate snakemake to get the correct python version**

# Annotate intersect using annovar by chromosome
```
$ qsub -W depend=afterok:17671637 /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/ANNOVAR/ANNOVAR_intersect_concat.pbs 
```

# Annotate concatenated intersect using SnpEff
```
$ qsub -W depend=afterok:17671637 /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/SnpEff/SnpEff_intersect_concat.pbs 
$ qsub -W depend=afterok:17671645 /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/SnpEff/SnpSift_filter_intersect_concat.pbs 
```

# Get union between annovar and snpeff
- Script to convert breed to the 10 target breeds and then other for genetic burden analysis.
```
$ convert_breeds.py
```
- Convert snpeff and annovar outputs to similar forms so that they can be compared
```
$ Get_annovar_snpeff_coding_tidy.py 
```
- Get union between annovar and snpeff
```
$ Get_union_annovar_snpeff_tidy.py
```
#Pull out exact intersect of coding variants for annovar/snpeff
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/bcftools_view.pbs 
```



