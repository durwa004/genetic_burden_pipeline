#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:11:06 2019

@author: durwa004
"""

import os
import gzip

#Get details about the type of variant,  etc.
            #Goals: Genes affected by GB
                #LOF variants
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/"
#path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/"

#Get list of genes to run through https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
#From Extract_lof_details.py
genes = {}
with open(path + "/lof_genes_with_symbols.txt", "r") as genes_file:
   genes_file.readline()
   for line1 in genes_file:
       line1 = line1.rstrip("\n").split("\t")
       a = line1[0].split(".")
       if "-" in line1[1]:
           if line1[0] in genes.keys():
               next
           else:
               genes[line1[0]] = line[0]
       else:
           if line1[0] in genes.keys():
               next
           else:
               genes[line1[0]] = line1[1]

#Add in gene constraint information
    HI = {}
    with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/disease_case_analysis/constraint_files/HI_Predictions_Version3.bed", "r") as input_file:
        input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            a = line[3].split("|")
            HI[a[0]] = line[4]

    DDG2P = {}
    with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/disease_case_analysis/constraint_files/DDG2P.csv", "r") as input_file:
        input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split(",")
            a = line[2].split('"')
            if len(a) == 1:
                DDG2P[line[0]] = a[0]
            elif len(a) == 3 or len(a) == 2:
                DDG2P[line[0]] = a[1]

    pLI = {}
    with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/disease_case_analysis/constraint_files/gnomad.v2.1.1.lof_metrics.by_transcript.txt", "r") as input_file:
        a = input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            if line[2] == "true":
                pLI[line[0]] = line[21]

#Get list of horse ids in order of vcf.
header = []
with gzip.open(path + "../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
            break

#Get breed info
horse_breed = {}
with open(path + "../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"

breeds = list(set(horse_breed.values()))

header = header[9:]

#Get list of variants from the annovar/snpeff intersect and just from annovar
variant_caller = {}
with open(path + "/lof_combined_intersect_lof_high_chrom_pos.txt", "r") as input_file:
#with open(path + "/snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[0] + ":" + line[1]
        variant_caller[ab] = line[2]

gb = {}
#with open(path + "/tabix_output/genetic_burden.txt", "r") as input_file, open(path + "/genetic_burden_tidy.txt", "w") as output_file:    
with open(path + "../gb_analysis/tabix_output/genetic_burden.txt", "r") as input_file, open(path + "/lof_tidy.txt", "w") as output_file, open(path + "/lof_tidy_brief.txt", "w") as output_b:
    print("VEP\tHI\tDDG2P\tpLI\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAC\tAF\tSNP\tconsequence\timpact\ttranscript\tgene\tcoding\tprotein", "\t".join(header), sep = "\t", file = output_file)
    print("VEP\tHI\tDDG2P\tpLI\tCHROM\tPOS\tREF\tALT\tAC\tAF\tSNP\tconsequence\timpact\ttranscript\tgene\tcoding\tprotein", sep = "\t", file = output_b)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        c_p = line[0] + ":" + line[1]
        if c_p in variant_caller.keys():
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
                consequence = "splice_region_variant"
            elif "5_prime" in consequence:
                consequence = "5_prime_UTR_variant"
            elif "frameshift" in consequence:
                consequence = "frameshift_variant"
            elif "stop_gained" in consequence:
                consequence = "stop_gained"
            elif "start_lost" in consequence:
                consequence = "start_lost"
            elif "stop_gained" in consequence:
                consequence = "stop_gained"
            elif "stop_lost" in consequence:
                consequence = "stop_lost"
            elif "gene_fusion" in consequence:
                consequence = "gene_fusion"
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
            transcript = gene
            transcript1 = transcript.split(".")
            transcript1 = transcript1[0]
            if transcript1 in genes.keys():
                gene=genes[transcript1]
            else:
                gene = transcript1
            if gene in HI.keys():
                HI_i = HI[gene]
            else:
                HI_i = "NA"
            if gene in DDG2P.keys():
                DDG2P_i = DDG2P[gene]
            else:
                DDG2P_i = "NA"
            if gene in pLI.keys():
                pLI_i = pLI[gene]
            else:
                pLI_i = "NA"
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
            print(variant_caller[c_p], HI_i, DDG2P_i, pLI_i,line[0], line[1], line[2], REF,ALT,line[5], line[6],AC, AF, SNP_se, consequence, impact, transcript, gene, coding, protein, "\t".join(line1), sep = "\t", file = output_file)
            print(variant_caller[c_p], HI_i, DDG2P_i, pLI_i,line[0], line[1], line[2], REF,ALT,line[5], line[6],AC, AF, SNP_se, consequence, impact, transcript, gene, coding, protein, sep = "\t", file = output_b)

#with open(path + "/lof_annovar.txt", "r") as annovar_file:#, open(path + "/lof_annovar_tidy.txt", "w") as output_file:
with open(path + "/genetic_burden_annovar.txt", "r") as annovar_file, open(path + "/genetic_burden_annovar_tidy.txt", "w") as output_file:
    print("VEP\tCHROM\tPOS\tREF\tALT\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein", "\t".join(header), sep = "\t", file = output_file)
    for line in annovar_file:
        line = line.rstrip("\n").split("\t")
        c_p = line[0] + ":" + line[1]
        if c_p in variant_caller.keys():
            REF = line[3]
            ALT = line[4]
            if "," in ALT:
                abcd = ALT.split(",")
                ALT = abcd[0]
            if len(REF) == len(ALT):
                SNP_ann = "SNP"
            else:
                SNP_ann = "indel"
            AF = line[8]
            if "," in AF:
                ef = AF.split(",")
                AF = ef[0]
            else:
                pass
            AC = line[7]
            if "," in AC:
                gh = AC.split(",")
                AC = gh[0]
            else:
                pass
            consequence = line[9]
            if "splice" in consequence:
                consequence = "splice_region_variant"
            elif "frameshift deletion" == consequence or "frameshift insertion" == consequence:
                consequence = "frameshift_variant"
            elif "stopgain" == consequence:
                consequence = "stop_gained"
            impact_ann = line[10]
            gene_ann = line[11].split("-")
            gene_ann = gene_ann[1]
            coding_ann = line[13]
            protein_ann = line[14].split(",")
            protein_ann = protein_ann[0]
            if line[10] == "LOW":
                next
            else:
                print(variant_caller[c_p], line[0], line[1], REF,ALT,AC, AF, SNP_ann, consequence, impact_ann,gene_ann, coding_ann, protein_ann, sep = "\t", file = output_file)

#Get details about shared variants
ann = {}
se = {}
with open(path + "/lof_annovar_tidy.txt", "r") as ann_in, open(path + "/lof_snpeff_tidy.txt", "r") as se_in, open(path + "/lof_details.txt", "w") as output_file, open(path + "/lof_details_brief.txt", "w") as output_b:
#with open(path + "/genetic_burden_annovar_tidy.txt", "r") as ann_in, open(path + "/genetic_burden_tidy.txt", "r") as se_in, open(path + "/genetic_burden_details.txt", "w") as output_file, open(path + "/genetic_burden_details_brief.txt", "w") as output_b:
    print("group\tchrom\tpos\tref\talt\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein\tSNP_ann\tconsequence_ann\timpact_ann\tgene_ann\tcoding_ann\tprotein_ann", "\t".join(header),sep = "\t",file = output_file)
    print("group\tchrom\tpos\tref\talt\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein\tSNP_ann\tconsequence_ann\timpact_ann\tgene_ann\tcoding_ann\tprotein_ann", file = output_b)
    ann_in.readline()
    se_in.readline()
    for line in ann_in:
        line = line.rstrip("\n").split("\t")
        a = line[1] + ":" + line[2]
        b = line[7] + ":" + line[8] + ":" + line[9] + ":" + line[10] + ":" + line[11] + ":" + line[12]
        if a in ann.keys():
            print(a)
        else:
            ann[a] = b
    for line in se_in:
        line = line.rstrip("\n").split("\t")
        ab = line[1] + ":" + line[2]
        if ab in se.keys():
            print(ab)     
        else:
            se[ab] = ":".join(line)
    variants = set.intersection(set(ann.keys()), set(se.keys()))
    variants = list(variants)
    for i,key in enumerate(variants):
         abc = se[key].split(":")
         abcd = ann[key].split(":")
         print("\t".join(abc[:13]), "\t".join(abcd), sep = "\t", file = output_b)
         print("\t".join(abc[:13]), "\t".join(abcd), "\t".join(abc[13:]), sep = "\t", file = output_file) 

#find variants unique to breeds  
unique = []
u_b = []
for item in gb.keys():
    a = gb[item].split(":")
    if len(set(a)) <3:
        unique.append(item)
        u_b.append(a[1])

with open(path + "/lof_tidy.txt", "r") as input_file, open(path + "/unique_lof.txt", "w") as output_file, open(path + "/unique_lof_brief.txt", "w") as output_b:
#with open(path + "/genetic_burden_details.txt", "r") as input_file, open(path + "/unique_gb.txt", "w") as output_file, open(path + "/unique_gb_brief.txt", "w") as output_b:
    header = input_file.readline()
    header = header.rstrip("\n").split("\t")
    print("breed", "\t".join(header), sep = "\t", file = output_file)
    print("breed", "\t".join(header[0:17]), sep = "\t", file = output_b)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[4] + ":" + line[5]
        for i in range(len(unique)):
            if a == unique[i]:
                if u_b[i] != "Other":
                    print(u_b[i], "\t".join(line), file = output_file, sep ="\t")
                    print(u_b[i],"\t".join(line[:20]),file = output_b, sep = "\t")
