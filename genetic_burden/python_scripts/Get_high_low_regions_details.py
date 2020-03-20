#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:11:06 2019

@author: durwa004
"""

import os
import gzip

#Get details about the type of variant,  etc.
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/"

#Get list of horse ids in order of vcf.
header = []
with gzip.open(path + "../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
            break

#Get breed info
horse_breed = {}
with open(path + "../../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"

breeds = list(set(horse_breed.values()))

header = header[9:]

#Need to get list of genes first
genes = {}
#with gzip.open(path + "Low_variation_regions_snpeff.coding.vcf.gz", "rt") as input_file:
with gzip.open(path + "Low_variation_regions_snpeff.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]: 
            next
        else:
            de = line[7].split("ANN=")
            bc = de[1].split("|")
            gene = bc[3]
            gene = gene.split("-")
            if "CHR_START" in gene[0]:
                gene = gene[-2]
            else:
                if "exon" not in gene[0]:
                    gene = gene[0]
                else:
                    if "id" in gene[1]:
                        gene = gene[2]
                    else:
                        gene = gene[1]
            genes[gene] = {}
#with gzip.open(path + "High_variation_regions_snpeff.coding.vcf.gz", "rt") as input_file:
with gzip.open(path + "High_variation_regions_snpeff.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            next
        else:
            de = line[7].split("ANN=")
            bc = de[1].split("|")
            gene = bc[3]
            gene = gene.split("-")
            if "CHR_START" in gene[0]:
                gene = gene[-2]
            else:
                if "exon" not in gene[0]:
                    gene = gene[0]
                else:
                    if "id" in gene[1]:
                        gene = gene[2]
                    else:
                        gene = gene[1]
            genes[gene] = {}

with open("High_low_variation_regions_all_genes.txt", "w") as output_file:
    for key in genes.keys():
        print(key, file = output_file)

#Add genes back in 
with open(path + "/bioDBnet_db2db_200131120203_835324622.txt", "r") as genes_file:
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       for key in genes.keys():
           b = key.split(".")
           if line1[1] == "-":
               if b[0] == a[0]:
                   genes[key] = a[0]
           elif b[0] == a[0]:
               genes[key] = line1[1]

gb = {}
#with gzip.open(path + "High_variation_regions_snpeff.coding.vcf.gz", "rt") as input_file, open(path + "/High_variation_regions_regions.txt", "w") as output_file,  open(path + "/High_variation_regions_regions_brief.txt", "w") as brief_f:
with gzip.open(path + "Low_variation_regions_snpeff.coding.vcf.gz", "rt") as input_file, open(path + "/Low_variation_regions_regions.txt", "w") as output_file,  open(path + "/Low_variation_regions_regions_brief.txt", "w") as brief_f:
#with gzip.open(path + "High_variation_regions_snpeff.vcf.gz", "rt") as input_file, open(path + "/High_variation_regions_all.txt", "w") as output_file,  open(path + "/High_variation_regions_all_brief.txt", "w") as brief_f:
#with gzip.open(path + "Low_variation_regions_snpeff.vcf.gz", "rt") as input_file, open(path + "/Low_variation_regions_all.txt", "w") as output_file,  open(path + "/Low_variation_regions_all_brief.txt", "w") as brief_f:
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein", sep = "\t", file = brief_f)  
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein", "\t".join(header), sep = "\t", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]: 
            next
        else:
            c_p = line[0] + ":" + line[1]
            ab = line[7].split(";")
            cd = ab[1].split("AF=")
            AF = cd[1]
            REF = line[3]
            ALT = line[4]
            if "," in ALT:
                abcd = ALT.split(",")
                ALT = abcd[0]
            if len(REF) == len(ALT):
                SNP_se = "SNP"
            else:
                SNP_se = "indel"
            if "," in AF:
                ef = AF.split(",")
                AF = ef[0]
            else:
                pass
            bc = ab[0].split("AC=")
            AC = bc[1]
            if "," in AC:
                gh = AC.split(",")
                AC = gh[0]
            else:
                pass
            gb[c_p] = AC
            de = line[7].split("ANN=")
            bc = de[1].split("|")
            consequence = bc[1]
            if "splice" in consequence:
                consequence = "splice_region"
            elif "5_prime" in consequence:
                consequence = "5_prime_UTR"
            elif "3_prime" in consequence:
                consequence = "3_prime_UTR"
            elif "downstream" in consequence:
                consequence = "downstream"
            elif "upstream" in consequence:
                consequence = "upstream"
            elif "intergenic" in consequence:
                consequence = "intergenic"
            elif "intragenic" in consequence:
                consequence = "intragenic"
            elif "intron" in consequence:
                consequence = "intronic"
            elif "non_coding_transcript" in consequence:
                consequence = "non_coding_transcript"
            elif "frameshift" in consequence:
                consequence = "frameshift"
            elif "stop_gained" in consequence:
                consequence = "stop_gain"
            elif "start_lost" in consequence:
                consequence = "start_loss"
            elif "stop_lost" in consequence:
                consequence = "stop_loss"
            elif "gene_fusion" in consequence:
                consequence = "gene_fusion"
            elif "conservative_inframe" in consequence or "disruptive_inframe" in consequence:
                consequence = "inframe_indel"
            elif "missense_variant" in consequence:
                consequence = "missense"
            elif "synonymous" in consequence:
                consequence = "synonymous"
            elif "initiator_codon_variant" in consequence:
                consequence = "initiator_codon"
            elif "stop_retained" in consequence:
                consequence = "stop_retained"
            coding = bc[9]
            protein = bc[10]
            impact = bc[2]
            gene = bc[3]
            gene = gene.split("-")
            if "CHR_START" in gene[0]:
                gene = gene[-1]
            else:
                if "exon" not in gene[0]:
                    gene = gene[0]
                else:
                    if "id" in gene[1]:
                        gene = gene[2]
                    else:
                        gene = gene[1]
            line1 = line[9:]
            for i in range(len(line1)):
                a = line1[i].split(":")
                if "/" in a[0]:
                    b = a[0].split("/")
                elif "|" in a[0]:
                    b = a[0].split("|")
                if b[0] == ".":
                    line1[i] = "missing"
                elif b[0] == "0":
                    if b[1] == "0":
                        line1[i] = "hom_WT"
                    if int(b[1]) > 0:
                        line1[i] = "het"
                        de = gb[c_p] + ":" + horse_breed[header[i]]
                        gb[c_p] = de
                elif int(b[0]) >0:
                    line1[i] = "hom"
                    de = gb[c_p] + ":" + horse_breed[header[i]]
                    gb[c_p] = de
            print(line[0], line[1], REF,ALT,AC, AF, SNP_se, consequence, impact, genes[gene], coding, protein, "\t".join(line1), sep = "\t", file = output_file)
            print(line[0], line[1], REF,ALT,AC, AF, SNP_se, consequence, impact, genes[gene], coding, protein, sep = "\t", file = brief_f)


#Get number of "LOC" genes in each type
high = 0
low = 0
count_h = 0
count_l = 0
with open(path + "/High_variation_regions_regions_brief_constraint.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        count_h +=1
        if "LOC" in line[0]:
            high +=1

with open(path + "/Low_variation_regions_regions_brief_constraint.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        count_l +=1
        if "LOC" in line[0]:
            low +=1
