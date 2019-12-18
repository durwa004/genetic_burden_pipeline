# genetic_burden_pipeline
Scripts and tools to estimate the genetic burden as part of my first aim of my thesis.
NB - search for **NEED TO WORK ON** for things that need to be done for different sections of thesis
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
$ python Extract_bcftools_by_ind_by_chrom.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools_without_Prze/ind_bcftools_stats_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -p bcftools
$ python Extract_bcftools_by_ind_by_chrom.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/ind_bcftools_stats_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -p gatk
```
Transfer to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/intersect_by_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_union/ind_bcftools_stats_files/union_by_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/*_AF_freq.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/AF_freq_files/
$ $ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools_without_Prze/ind_bcftools_stats_files/bcftools_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/ind_bcftools_stats_files/gatk_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
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
$ GB_by_individual.R
```
#Get exact breed details for looking at GB variants
```
$ Get_genetic_burden_by_individual_breed.py
```
#Get additional details from gb e.g. genes involved etc.
```
$ Get_genetic_burden_details.py
```
#Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/ann_se_gb_by_individual.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/ann_se_lof_by_individual.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/genetic_burden_535_horses.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/genetic_burden_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/lof_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/unique_gb.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ GB_by_individual.R
$ Extract_GB_details.py
```

# Need to find possible DCVs (probably falls in lof analysis)  
- Pull out all variants with no homozygotes
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_no_homozygous_sites.pbs/ 
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_variant_details.py
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/no_homozygotes/no_homozygotes_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/no_homozygotes/lof_variants_no_homozygotes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/no_homozygotes/unique_no_homozygotes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/no_homozygotes/thesis_intersect_no_homozygotes.vcf.gz.stats /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
```
# Get variants present in all individuals
- All homozygotes
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_all_homozygous_sites.pbs 
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/view_no_homozygous_sites.pbs/Get_variant_details.py
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/all_homozygotes_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/lof_variants_all_homozygotes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/thesis_intersect_all_homozygotes.vcf.gz.stats /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
```
- All heterozygous/homozygous
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_af_over_0.5.pbs 
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_variants_present_in_all.py
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_intersect_variants_present_in_all.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_intersect_variants_present_in_all.vcf.gz 
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/variants_present_in_all_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/variants_present_in_all_not_homozygous_in_all.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/lof_variants_present_in_all.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/thesis_intersect_variants_present_in_all.vcf.gz.stats /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
```
- Analyse data
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/Extract_homozygous_heterozygous_details.py
$ homozygous_heterozygous_analysis.R
```
# Compare allele frequency of gb variants and non-gb variants
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_AF_all_variants.py
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/AF_gb_cf_all_variants/AF_gb_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/AF_gb_cf_all_variants/AF_all_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ Compare_AF_gb_all_variants.R
```

# Number of variants unique to populations
Plan to split the thesis intersect by breed group
# Number of variants with big differences in frequency between populations (<3% in one and >10% in another)
- Get breed vcfs
```
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_shell.sh
```
- Get common/rare variants from each breed
```
$ qsub -I -l nodes=1:ppn=8,walltime=12:00:00,mem=4g
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_rare_common_af_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/ -r 0.03 -c 0.10

$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {breed}_rare.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {breed}_rare.vcf.gz 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {breed}_common.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {breed}_common.vcf.gz 
```
- Get intersect between the files
```
$ python ../../python_scripts/Generate_bcftools_isec_rare_common_variants_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/ -b WP
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_isec_extract_rare_common_variants/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_isec_extract_rare_common_variants/pbs_shell.sh
```
- Creates a lot of spare files and lots of GB! So tidy up files and get bcftools stats information from them
```
$ delete_unnecessary_files.py 
$ sh breed_common_rare_variants_tidy.sh 
```
- Get chrom/pos list (no header) to extract variants that are rare in breed and common in other breed (0002.vcf.gz) these variants from the SnpEff file
```
$ python ../../python_scripts/Extract_bcftools_stats_isec.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/
$ sed -i 's$ xxxxx$$g' *
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_query_breed_rare_common/bcftools_query_chrom_pos.sh 
```
- Extract the breed_rare_common variants from the snpeff file
```
$ python ../../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_chrom_pos/ 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions/Extract_variants_rare_common_breed.pbs
$ module load bcftools
$ bcftools stats rare_common_breed_snpeff.vcf.gz > rare_common_breed_snpeff.vcf.gz.stats
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_breed_rare_common.pbs
$ cat rare_common_breed_high.vcf.gz rare_common_breed_moderate.vcf.gz rare_common_breed_low.vcf.gz > rare_common_breed_coding.vcf.gz
```
- Extract the breed details from the output vcf for breed rare/common variants
- Need to pull in info re which breed rare/common
```
$ qsub -I -l nodes=1:ppn=8,walltime=12:00:00,mem=40g
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_rare_common_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_snpeff -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_chrom_pos/
```
- Move files back to my laptop 
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_snpeff/breed_rare_other_breed_common.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_snpeff/breed_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ Get_breed_rare_breed_common_details.py
```
- Get genes details
https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
```
$ split -l 2,500 breed_breed_disc_gene_list.txt breed_breed_genes_
```
https://david.ncifcrf.gov/gene2gene.jsp

# Number of variants that are rare in one breed and common in the general population and variants that are common in that breed and rare in the general population 

- Create vcfs without each breed
```
$ python ../../python_scripts/Generate_bcftools_view_extract_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_view_remove_breed/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_remove_breeds/pbs_shell.sh
```
- Create rare AF <3% and common AF >10% pop vcfs
```
$ qsub -I -l nodes=1:ppn=8,walltime=24:00:00,mem=4g -q lab-long
$ source activate snakemake
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_rare_common_breed_pop_af.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/ -r 0.03 -c 0.10

$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip without_{breed}_rare.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix without_{breed}_rare.vcf.gz 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip without_{breed}_common.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix without_{breed}_common.vcf.gz 
```
- Get intersect between the files
```
$ python ../../python_scripts/Generate_bcftools_isec_rare_common_variants_breed_pop.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt 
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_isec_breed_pop_rare_common_variants/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_isec_breed_pop_rare_common_variants/pbs_shell.sh
```
- Creates a lot of spare files and lots of GB! So tidy up files 
```
$ delete_unnecessary_files.py 
$ module load bcftools
$ sh breed_pop_common_rare_tidy.sh 
```
- Get chrom/pos list (no header) to extract variants that are rare in breed and common in other breed (0002.vcf.gz) these variants from the SnpEff file
```
$ python ../../python_scripts//Extract_bcftools_stats_isec_po.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/
$ sed -i 's$ xxxxx$$g' *
$ qsub -I -l nodes=1:ppn=1,walltime=12:00:00,mem=4g
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_query_breed_rare_common/bcftools_query_pop_chrom_pos.sh 
```
- Extract the breed_rare_common variants from the snpeff file
```
$ python ../../python_scripts/Generate_extract_variants_pop.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_breed_common_rare_pop.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_breed_rare_common_pop.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_unique.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_breed_rare_common_pop.pbs
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_breed_common_rare_pop.pbs
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_unique.pbs
```
- Extract the breed details from the output vcfs for breed rare/common variants **LOOK AT NEXT**
```
$ qsub -I -l nodes=1:ppn=1,walltime=12:00:00,mem=4g
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_pop_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ -b unique
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_pop_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ -b rare
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_pop_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ -b common
```
- Transfer to my laptop for analysis
```
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_unique_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_unique_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_rare_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_rare_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_common_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_common_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/Get_breed_pop_rare_common_details.py
```
- Get genes details
https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
https://david.ncifcrf.gov/gene2gene.jsp



# known causal variants 
Decided to do all variants on OMIA - download from: https://omia.org/results/?search_type=advanced&gb_species_id=9796 (12/11/2019) - 154 variants (also PSSM2 patent, and 88 horses paper) - IDENTIFIED 93 OF THEM
- Also pulled equine qtls from animal genome: https://www.animalgenome.org/cgi-bin/QTLdb/EC/download?tmpname=mapDwnLd&file=cM
As of Release 39, there have been 2,260 horse QTLs released for public access on the Horse QTLdb. These data were curated from 88 publications and represent 54 different horse traits.
Excluded QTLs without snp IDs/positions on the genome in the download
Selected the peak SNP to investigate = 1,730 QTLs (1,559 unique)
679 mapped uniquely, 549 mapped to multiple places in the genome
633 were present in at least one horse
- Tidy up variant locations (from my computer and create shell script)
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Tidy_known_variants_for_extraction.py 
```
- Get variant locations - on my laptop
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py
```
- Pull out QTLs from dbsnp **WORKING ON**
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Tidy_known_variants_for_extraction.py
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/grep_get_rs_chrom_pos.sh durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/known_qtls/
$ sh grep_get_rs_chrom_pos.sh
$ source activate ensembl-vep
$ split -l 500 QTL_EC2_chrom_pos.bed QTLs_
$ sh ncbi_remap_known.sh
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants_pop.py
$ sh QTLs_remapped.sh
$ qsub -I -l nodes=1:ppn=1,mem=2g,walltime=24:00:00
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants_pop.py
```
- Move back to my laptop 
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_QTLs_from_tabix/known_QTL_locations/known_QTLs_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_QTLs_from_tabix/known_QTL_locations/No_QTLs_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/Get_correct_genotypes_for_known_variants.py
$ known_causal_variants.R
```

- Extract known variants - MSI
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/known_SNVs_tabix.sh durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/known_SNVs.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/
$ sh known_SNVs_tabix.sh
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py 
```
- Move back to my laptop - Manually check position of variants (and get combination genotypes for the variants that a group of variants represent the phenotype
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations/known_variants_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations/No_variants_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/Get_correct_genotypes_for_known_variants.py
$ known_causal_variants.R
```

# Future things to look at
--Fst of genes containing these variants - look for genes with strong differentiation between populations
