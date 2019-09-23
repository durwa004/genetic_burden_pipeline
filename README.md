# genetic_burden_pipeline
Scripts and tools to estimate the genetic burden as part of my first aim of my thesis.

**Input = Intersect of bcftools/gatk haplotype caller (group calling) on all NC_ chromosomes (exclude unplaced for now).**

NB - will need to re-estimate for the individual calls as well.

**All python scripts that run from command line need source activate snakemake to get the correct python version**

# Annotate intersect using annovar
```
$ python ../python_generation_scripts/Generate_ANNOVAR_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/thesis_workflow/VARIANT_ANNOTATION/ANNOVAR/ANNOVAR_intersect_NC_009149_3.pbs 
```
**Concatenate annovar output files**
```
$ python ../python_generation_scripts/Generate_cat_annovar.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/annovar/ 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/ANNOVAR/annovar_cat.pbs
```

# Annotate intersect using SnpEff
```
$ python ../python_generation_scripts/Generate_SnpEff_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/thesis_workflow/VARIANT_ANNOTATION/SnpEff/SnpEff_intersect_NC_009149_3.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/thesis_workflow/VARIANT_ANNOTATION/SnpEff/SnpSift_filter_NC_009149_3.pbs 
```
**Concatenate SnpEff output files**
```
$ python ../python_generation_scripts/Generate_bcftools_concat.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/SnpEff
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/thesis_workflow/VARIANT_ANNOTATION/SnpEff/snpeff_concat.pbs 
```

# Get union between annovar and snpeff
**Get chrom/pos/impact of all coding variants**
- Also creates info file with details on numbers of variants
```
$ python Get_chrom_pos_impact.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -p SnpEff
$ python Get_chrom_pos_impact.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -p annovar
```
**Output chrom/pos of union and intersect of the annotation programs**
- Also creates an info file with details on numbers of variants
```
$ python Get_union_intersect_VEPs.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/
```

# Miscellanous scripts
- Script to convert breed to the 10 target breeds and then other for genetic burden analysis.
```
$ convert_breeds.py
```
