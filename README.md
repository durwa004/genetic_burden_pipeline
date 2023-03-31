# genetic burden pipeline
Scripts and tools to estimate the genetic burden as part of my first aim of my thesis.

#Input vcfs (SNPs and indels): /panfs/jay/groups/6/durwa004/shared/PopulationVCF

#Remove two individuals with very high GB (>3 SDs higher than the mean)

```
$ gatk SelectVariants -V /home/durwa004/shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz --exclude-sample-name M989 --exclude-sample-name M6468 -O joint_genotype_subset.vep.vcf
```

# Get tidy breeds (i.e. with other)
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/convert_breeds.py
```

# Get the number of variants, tstv, and DOC for each horse
```
$ sbatch ../scripts/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_stats_per_horse.slurm
$ python ../scripts/genetic_burden_pipeline/genetic_burden/python_scripts/Extract_bcftools_stats_ind.py -d bcftools_stats_combined_per_horse.stats
```

##Move DOC and no. reads back to my laptop
```
##$ scp durwa004@mesabi.msi.umn.edu:/home/durwa004/durwa004/genetic_burden/ind_number_of_variants.txt
##$ bcftools_stats_analysis.R
```

#Get number of variants
```
$ bcftools stats ../../shared/PopulationVCF/joint_genotype_indels.goldenPath.vep.vcf.gz > indels.stats
$ bcftools stats ../../shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz > SNPs.stats
```

#Run SnpEff
```
$sbatch /home/durwa004/durwa004/scripts/genetic_burden_pipeline/variant_annotation/SnpEff/SnpEff.slurm 
```

#Pull out high/moderate/low impact variants
```
$ sbatch /home/durwa004/durwa004/scripts/genetic_burden_pipeline/variant_annotation/SnpEff/SnpSift.slurm 
$ sbatch /home/durwa004/durwa004/scripts/genetic_burden_pipeline/variant_annotation/Ensembl-VEP/Ensembl-VEP_filter.pbs
```

#Pull out type of variant
```
$ python /home/durwa004/durwa004/scripts/genetic_burden_pipeline/genetic_burden/python_scripts/Get_type_of_variant.py -d joint_genotype_combined.goldenPath.snpeff.hml.vcf.gz -p SnpEff
$ python /home/durwa004/durwa004/scripts/genetic_burden_pipeline/genetic_burden/python_scripts/Get_type_of_variant.py -d joint_genotype_combined.goldenPath.vep.hml.vcf.gz -p VEP
```

#Get intersect between VEP and SnpEff
```
$  /home/durwa004/durwa004/scripts/genetic_burden_pipeline/genetic_burden/python_scripts/Get_intersect_SnpEff_VEP.py
```

#Transfer to my laptop and analyze
```
$ scp durwa004@mesabi.msi.umn.edu:/home/durwa004/durwa004/genetic_burden/SnpEff_VEP_intersect.txt GB_project
$ scp durwa004@mesabi.msi.umn.edu:/home/durwa004/durwa004/genetic_burden/SnpEff.VEP.intersect.txt GB_project
$ /home/durwa004/durwa004/scripts/genetic_burden_pipeline/R_analysis/GB_paper.R
```

#Get number of variants by individual
```
$ Get_genetic_burden_by_individual.py
$ scp durwa004@mesabi.msi.umn.edu:/home/durwa004/durwa004/genetic_burden/SnpEffVEP.intersect.individual.txt GB_project
$  /home/durwa004/durwa004/scripts/genetic_burden_pipeline/R_analysis/GB_paper.R
```

#Outliers - thesis union
```
$ Extract_bcftools_stats_ind.py
$ gatk SelectVariants -R /home/durwa004/durwa004/GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fasta -V thesis_union.vcf.gz -O thesis_union.intersect.vcf -L /home/durwa004/durwa004/genetic_burden/SnpEff.VEP.in
tersect.pos.list    
$ Get_genetic_burden_by_ind.py
$ scp durwa004@mesabi.msi.umn.edu:/scratch.global/marlo072/CheckHorses/thesis_union/tu.SnpEff.VEP.intersect.individual.txt GB_project
```
#Outliers - thesis intersect
```
$ s3cmd sync s3://durwa004_2019/thesis_analysis/thesis_intersect/thesis_intersect/thesis_intersect.vcf.gz ../thesis_intersect/                                                                                     
$ Extract_bcftools_stats_ind.py
$ gatk SelectVariants -R /home/durwa004/durwa004/GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fasta -V thesis_intersect.vcf.gz -O thesis_intersect.intersect.vcf.gz -L /home/durwa004/durwa004/genetic_burden/SnpEff.VEP.intersect.pos.list
$ Get_genetic_burden_by_ind.py
$ scp durwa004@mesabi.msi.umn.edu:/scratch.global/marlo072/CheckHorses/thesis_intersect/ti.SnpEff.VEP.intersect.individual.txt GB_project
```

####Remove M989 and M6468
```
$ sbatch ../scripts/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_subset.slurm 
```
#Then run through steps above

###OLD###
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
$ python ../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/tabix_output -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/thesis_intersect_snpeff.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/Extract_variants.pbs 
$ cat *.txt > genetic_burden.txt
$ python ../python_scripts/Generate_extract_annovar_variants.py  -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/annovar/annovar_coding_tidy.txt -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt 
```
#Get genetic burden
```
$ Get_genetic_burden_details.py
```
#Figure out genetic burden per individual/breed
```
$ Get_genetic_burden_by_individual.py
$ GB_by_individual.R
```
Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_by_individual.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
```
#Get exact breed details for looking at GB variants e.g. genes involved etc.
#Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_details_brief.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/unique_gb_brief.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper
$ sed -i 's/#/_/g' genetic_burden_details_brief.txt
$ GB_by_individual.R
$ Extract_GB_details.py
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_genes.txt /Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ split -l 500 genetic_burden_genes.txt gb_genes_
```
#Run genes on https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
```
$ cat /Users/durwa004/Downloads/bioDBnet_db2db_20071610* > genetic_burden_genes_with_symbols.txt
$ scp /Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/genetic_burden_genes_with_symbols.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/
$ Extract_GB_details.py
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_genes_AF_over_50.txt /Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_genes_details.txt /Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_genes_over_5_variants.txt /Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_details_brief_tidy.txt /Users/durwa004/Documents/Postdoc/Phd_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ GB_gene_analysis.R
```

- Pull out just lof variants
Get union between annovar and snpeff
```
$ Get_union_annovar_snpeff_tidy.py
```
Extract variants from annovar and snpeff
```
$ python ../python_scripts/Generate_extract_snpeff_variants.py  -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/thesis_intersect_snpeff.coding.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_combined_intersect_lof_high_chrom_pos.txt
$ python ../python_scripts/Generate_extract_annovar_variants.py  -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/annovar/annovar_coding_tidy_lof.txt -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_combined_intersect_lof_high_chrom_pos.txt
```
#Figure out LOF per individual/breed
```
$ Get_genetic_burden_by_individual.py
```
Move to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_by_individual.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/nature_genetics_paper/
```
Analyse
```
$ lof_by_individual.R
```
#Get additional details from gb e.g. genes involved etc.
```
$ Get_genetic_burden_details.py
$ Extract_lof_details.py
```
#Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_details.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_details_brief.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_details.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/genetic_burden_details_brief.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/unique_lof.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/unique_lof_brief.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/gb_analysis/
$ sed -i -- 's/#/_/g' lof_details*
$ lof_by_individual.R
```
Look at gene info
```
$ Extract_lof_details.py
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_genes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/
```
Get gene details: #Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol
Move back to MSI for analysis
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/lof_genes_with_symbols.txt  durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/
$ Extract_lof_details.py
```
Move high AF genes and genes with multiple LOF variants to my laptop for gene enrichment
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_multi-variant_genes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/lof_high_AF_genes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/lof/
```
DAVID 6.8 (use horse and human as species) - results in poster = functional clustering
https://david.ncifcrf.gov/summary.jsp

#Need to convert gene symbols to Gene description to convert LOC1 genes
Get gene details: #Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #gene symbol to Gene description
        #Helpful info https://www.biostars.org/p/345074/
```
$ Get_lof_gene_unfo.py
$ sh lof_high_AF_for_Uniprot.txt
$ Get_lof_gene_unfo.py
```
        http://www.ensembl.org/biomart/martview
        Details: https://www.biostars.org/p/179914/
        
Get number of common and rare variants GB and LOF variants.
Then take 1% of each to determine % false positives

# Validate 5% of GB variants (make sure it is 5% GB only and 5% LOF/GB variants)
```
$ genetic_burden/python_scripts/Get_variants_for_false_postive_analysis.py 
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


# known causal variants #
Decided to do all variants on OMIA - download from: https://omia.org/results/?search_type=advanced&gb_species_id=9796 
05/14/2020 - 46 phenotypes have causative variants reported
  - Download this table (copy and paste), then order by gene name and click on each OMIA ID and download each .csv file individually.
  - 239 genetic traits and disorders
  - 97 variants reported as causative for these 46 phenotypes.
  - 14 variants (11 phenotypes) were reported as gross (>20 bp) structural variants and were excluded
  - 83 variants (38 phenotypes) were SNPs or small (<= 20 bp) structural variants for investigation.
  - 29/83 were present in this dataset
  
NB - add in PSSM2 patent, and 88 horses paper

```
$ scp /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/OMIA_variants/OMIA_SNVs_05_14_20.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_variants_May_2020/
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Tidy_known_variants_for_extraction.py 
```
- Get variant locations
```
$ sh ../known_SNVs_tabix.sh
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py
```
- Move back to my laptop to double check locations
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_variants_May_2020/No_variants_present.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/OMIA_variants/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_variants_May_2020/known_variants_present.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/OMIA_variants/
```
- Move back to MSI for additional analysis
```
$ scp /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/OMIA_variants/known_variants_present_exact_locations_May_2020.txt  durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_variants_May_2020/
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py
```
- Move back to my laptop for analysis
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/variants_by_individual_R.txt /Users/durwa004/Documents/Postdoc/PhD_papers_for_publication/Nature_genetics/Post_thesis/OMIA_variants/
$ known_causal_variants.R
```

- Also pulled equine qtls from animal genome: https://www.animalgenome.org/cgi-bin/QTLdb/EC/download?tmpname=mapDwnLd&file=cM
As of Release 39, there have been 2,260 horse QTLs released for public access on the Horse QTLdb. These data were curated from 88 publications and represent 54 different horse traits.
Excluded QTLs without snp IDs/positions on the genome in the download
Selected the peak SNP to investigate = 1,730 QTLs (1,559 unique)
679 mapped uniquely, 549 mapped to multiple places in the genome
633 were present in at least one horse
- Tidy up variant locations (from my computer and create shell script)
- For some reason the OMIA download file got screwed up and so created new table to rerun this (known_variants_tables.xls)
Move variant locations to MSI
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/known_variants_locations.txt  durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_locations_2020/
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Tidy_known_variants_for_extraction.py 
```
- Get variant locations
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py
```
- Move back to my laptop to double check locations
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/No_variants_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/known_variants_present.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/
```
- Send back to MSI to get number per indidivual
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/known_variants_present_exact_locations.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/
```
- Move back to my laptop for analysis
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/variants_bt_indvidual.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/
$ known_causal_variants.R
```
- Get tidy table of known disease causing variants - Mick wants me to stop using CC_ so change to coat_col_
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Tidy_known_variants_for_extraction.py 
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/known_variants_present_exact_locations_tidy_april_2020.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/
```
- Convert by individual table into R friendly version 
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_known_causal_variants.py
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_disease_locations_2020/variants_by_indvidual_R.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/
```

- Pull out QTLs from dbsnp
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
#Get final variant details (from table that is in paper)
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/animalgenomeQTL_for_extraction.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/
$ Get_known_causal_variants_pop.py
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/QTLs_table_tidy.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/2020/
```

