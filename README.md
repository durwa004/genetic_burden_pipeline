# genetic_burden_pipeline
Scripts and tools to estimate the genetic burden as part of my first aim of my thesis.

To do for paper: 
- Each variant caller, across population, breed, variant location (exon, intron, intergenic, splice site, or within 5kb of a gene), SNP consequence (high/moderate/low/modifier).
  - TsTv
  - HetNRhom ratio
  
# Get tidy breeds (i.e. with other)
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/convert_breeds.py
```

# Get DOC and number of reads
```
$ python python_scripts/Generate_unfreeze_coverage.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/doc/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -c get -ft coverage.tsv
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/s3_scripts/s3cmd_sync_coverage.tsv.pbs 
$ python Extract_coverage.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/doc/coverage/
```
Move DOC and no. reads back to my laptop
```
$ scp durwa004@login02.msi.umn.edu://home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/doc/coverage/summary_files/* ../DOC
$ bcftools_stats_analysis.R
```

# Number of variants called ST/GATK and intersect/union
**bcftools_stats**
- Run bcftools for all files in a directory
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/Generate_bcftools_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -e _intersect.vcf.gz
```
- Extract bcftools info
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/Extract_bcftools_stats.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -p bcftools
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_stats.genotyped.vcf.gz.pbs 
```
- Now I have the concatenated intersect
```
$ bcftools stats thesis_intersect.vcf.gz > thesis_intersect.vcf.gz.stats
```
- Move stats info back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/thesis_intersect.vcf.gz.stats /Users/durwa004/Desktop
```
- Analyze stats and get tstv ratio etc
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/R_analysis/bcftools_stats_analysis.R
```

# Determine number of variants per individual - output file should also have breed information
- Need to do for union/intersect
```
$ python ../Generate_bcftools_by_individual.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -ind /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt 
$ python ../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_by_ind/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_by_ind/pbs_shell.sh 
$ python Extract_bcftools_by_ind_stats.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/
```
- For gatk/bcftools
```
$ python ../../python_scripts/Generate_bcftools_by_chr_by_ind.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/genotyped_files/bcftools/ -e .genotyped.vcf.gz -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt 
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/pbs_scripts/ind_gatk_bcftools/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/pbs_scripts/ind_gatk_bcftools/pbs_shell.sh 
```
NEED TO WORK ON
```
$ python ../python_scripts/Extract_bcftools_by_ind_by_chrom.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/
```
Transfer to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/intersect_by_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/*_AF_freq.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/AF_freq_files/
```
-- Breed differences in the number of variants per individual
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/R_analysis/bcftools_stats_analysis.R
```
-- AF of observed variants
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/merge_bcfstats_AF_191015.py
```

**Input = Intersect of bcftools/gatk haplotype caller (group calling) on all chromosomes**

NB - will need to re-estimate for the individual calls as well.

**All python scripts that run from command line need source activate snakemake to get the correct python version**

# SnpEff/Annovar analysis
#Annotate intersect using annovar by chromosome (won't work on concatenated file)
#If need to get by chromosome from intersect file: qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/ANNOVAR/bcftools_view_by_chr.pbs 

```
$ python ../python_generation_scripts/Generate_ANNOVAR_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/
$ python ../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/ANNOVAR/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/ANNOVAR/pbs_shell.sh 
$ cat * >thesis_intersect.exonic_variant_function
$ cat * >thesis_intersect.variant_function

```

#Annotate concatenated intersect using SnpEff
```
$ qsub -W depend=afterok:17671637 /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/SnpEff/SnpEff_intersect_concat.pbs 
$ qsub -W depend=afterok:17671645 /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/variant_annotation/SnpEff/SnpSift_filter_intersect_concat.pbs 
```

#Pull out type of variant - all
```
$ Get_type_of_variant_SnpEff.py 
$ Get_type_of_variant_annovar.py 
```
Transfer to my laptop and analyze
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/annovar/annovar_variant_type_all.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/variant_type_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/SnpEff_variant_type_all.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/variant_type_analysis/
$ variant_type_analysis.R
```

#Pull out type of variant - coding
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_type_of_variant_SnpEff_coding.py 
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_type_of_variant_annovar_coding.py
```
Transfer to my laptop and analyze
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/SnpEff_variant_type_coding.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/variant_type_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/annovar/annovar_variant_type_coding.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/variant_type_analysis/
$ variant_type_analysis.R
```
#Get union between annovar and snpeff
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
#Pull out combined intersect of high/moderate variants for annovar/snpeff
```
$ python ../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/thesis_intersect_snpeff.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt 
$ qsub /home/mccuem/
shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/Extract_variants.pbs 
$ cat HIGH_HIGH* > ann_se_high_high.txt
$ cat HIGH_MODERATE* > ann_se_high_moderate.txt
$ cat MODERATE_HIGH* > ann_se_moderate_high.txt
$ cat ann_se_* > genetic_burden.txt
```
#Figure out genetic burden per individual/breed
```
$ Get_genetic_burden_by_individual.py
```
#Get additional details from gb e.g. genes involved etc.
```
$ Get_genetic_burden_details.py
```
#Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/ann_se_gb_by_individual.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/genetic_burden_535_horses.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/genetic_burden_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/lof_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/unique_gb.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ GB_by_individual.R
$ Extract_GB_details.py
```

# Need to find possible DCVs (probably falls in lof analysis)  
# Number of variants unique to populations
Plan to split the thesis intersect by breed group
- Number of variants with big differences in frequency between populations (<3% in one and >10% in another)
```
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b Shetland -m 3 -c 11
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b TB -m 3 -c 10
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b STB -m 3 -c 10
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b Icelandic -m 1 -c 3
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b WP -m 1 -c 4
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b Morgan -m 1 -c 4
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b Belgian -m 1 -c 4
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b Clydesdale -m 1 -c 4
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b QH -m 4 -c 14
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -b Arabian -m 2 -c 7
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_view_extract_rare_common_variants/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_shell.sh
$ python ../../python_scripts/Generate_bcftools_isec_for_rare_common_variants_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -b Arabian
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_isec_extract_rare_common_variants/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_isec_extract_rare_common_variants/pbs_shell.sh
```
Creates a lot of spare files and lots of GB! So tidy up files
```
$ delete_unnecessary_files.py 
$ sh breed_common_rare_variants_tidy.sh 
```

- Number of variants that are rare in one breed and common in the general population and variants that are common in that breed and rare in the general population 
#Tried AF cut off for rare population of 0.05% and common population 5%
#Need to figure out how to run this to exclude the horses of that breed in the population AF
```
$ python ../../Generate_bcftools_view_rare_variants_breed_common_variants_pop.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -b Arabian -m5 4
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_rare_variants_breed_common_variants_pop/
$  sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_rare_variants_breed_common_variants_pop/pbs_shell.sh
```
Then need to pull out the population allele frequencies using GATK
```
$ python ../../python_scripts/Generate_gatk_extract_by_pop_AF.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_rare_variants_breed_common_variants_pop/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_rare_variants_breed_common_variants_pop/pbs_shell.sh
$ Get_rare_breed_common_pop_info.py 
```
Transfer to my laptop and analyze using R, then extract list of chrom/pos (unique) to extract from SnpEff file
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/rare_breed_common_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/common_breed_rare_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ breed_pop_variants.R
$ Theses_analysis.py
```
Then transfer list of chrom/pos for the variants to be extracted from SnpEff output
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/rare_breed_common_pop_shared_variants.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/common_breed_rare_pop_shared_variants.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/
```
Extract information about breed/pop difference variant type
```
$ python ../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/rare_breed_common_pop -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/SnpEff/thesis_intersect_snpeff.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/rare_breed_common_pop_shared_variants.txt 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/Extract_variants.pbs 
$ cat * > rare_breed_common_pop_snpeff_variants.txt
$ python ../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/common_breed_rare_pop -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/SnpEff/thesis_intersect_snpeff.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/common_breed_rare_pop_shared_variants.txt 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/Extract_variants.pbs
$ cat * > common_breed_rare_pop_snpeff_variants.txt
$ cat * > rare_breed_common_pop_snpeff_variants.txt
$ Extract_breed_pop_shared_variant_information.py
```
Download back to my laptop for analysis!
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/common_breed_rare_pop/common_breed_rare_pop_snpeff_info.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/rare_breed_common_pop/rare_breed_common_pop_snpeff_info.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
```
--Fst of genes containing these variants - look for genes with strong differentiation between populations
-Putatively functional variants (high/moderate) - number per individual genome
-- breed differences
- Predicted consequence of variant

# known causal variants
- Tidy up variant locations (from my computer and create shell script)
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Tidy_known_variants_for_extraction.py 
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/known_SNVs_tabix.sh
```
- Get variant locations - on my laptop
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py
```
- Extract known variants - MSI
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/known_SNVs_tabix.sh durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/known_SNVs.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/
$ sh known_SNVs_tabix.sh
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py 
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations/known_variants_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants
$ known_causal_variants.R
```


