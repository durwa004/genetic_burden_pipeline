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
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[0] + ":" + line[1]
        variant_caller[ab] = line[2]

gb = {}
with open(path + "/lof_snpeff.txt", "r") as input_file, open(path + "/lof_snpeff_tidy.txt", "w") as output_file:    
    print("VEP\tCHROM\tPOS\tREF\tALT\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein", "\t".join(header), sep = "\t", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
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
        for i in range(len(line)):
            if "0/1" in line[i] or "0/2" in line[i] or "0/3" in line[i]:
                line[i] = "het"
            elif "1/1" in line[i] or "2/2" in line[i] or "1/2" in line[i] or "1/3" in line[i] or "2/3" in line[i] or "2/1" in line[i]:
                line[i] = "hom"
            elif "./." in line[i]:
                line[i] = "missing"
            elif "0/0" in line[i]:
                line[i] = "ref"
        print(variant_caller[c_p], line[0], line[1], REF,ALT,AC, AF, SNP_se, consequence, impact, gene, coding, protein, "\t".join(line[9:]), sep = "\t", file = output_file)
        line1 = line[9:]
        for i in range(len(line1)):
            if "het" in line1[i] or "hom" in line1[i]:
                de = gb[c_p] + ":" + horse_breed[header[i]]
                gb[c_p] = de

with open(path + "/lof_annovar.txt", "r") as annovar_file, open(path + "/lof_annovar_tidy.txt", "w") as output_file:
    print("VEP\tCHROM\tPOS\tREF\tALT\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein", "\t".join(header), sep = "\t", file = output_file)
    for line in annovar_file:
        line = line.rstrip("\n").split("\t")
        c_p = line[0] + ":" + line[1]
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
        for i in range(len(line)):
            if "0/1" in line[i] or "0/2" in line[i] or "0/3" in line[i]:
                line[i] = "het"
            elif "1/1" in line[i] or "2/2" in line[i] or "1/2" in line[i] or "1/3" in line[i] or "2/3" in line[i] or "2/1" in line[i]:
                line[i] = "hom"
            elif "./." in line[i]:
                line[i] = "missing"
            elif "0/0" in line[i]:
                line[i] = "ref"
        print(variant_caller[c_p], line[0], line[1], REF,ALT,AC, AF, SNP_ann, consequence, impact_ann,gene_ann, coding_ann, protein_ann,"\t".join(line[9:]), sep = "\t", file = output_file)

#Get details about shared variants
ann = {}
se = {}
ij = 0
with open(path + "/lof_annovar_tidy.txt", "r") as ann_in, open(path + "/../annovar/annovar_coding_tidy.txt", "r") as ann_tidy, open(path + "/lof_snpeff_tidy.txt", "r") as se_in, open(path + "/lof_details.txt", "w") as output_file, open(path + "/lof_details_brief.txt", "w") as output_b:
    print("group\tchrom\tpos\tref\talt\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein\tSNP_ann\tconsequence_ann\timpact_ann\tgene_ann\tcoding_ann\tprotein_ann", "\t".join(header),sep = "\t",file = output_file)
    print("group\tchrom\tpos\tref\talt\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein\tSNP_ann\tconsequence_ann\timpact_ann\tgene_ann\tcoding_ann\tprotein_ann", file = output_b)
    ann_in.readline()
    se_in.readline()
    for line in ann_in:
        line = line.rstrip("\n").split("\t")
        a = line[1] + ":" + line[2]
        b = line[7] + ":" + line[8] + ":" + line[9] + ":" + line[10] + ":" + line[11] + ":" + line[12]
        ann[a] = b
    for line in se_in:
        line = line.rstrip("\n").split("\t")
        ab = line[1] + ":" + line[2]
        if ab in ann.keys():
            cd = ann[ab].split(":")
            print("\t".join(line[:13]), "\t".join(cd), "\t".join(line[13:]), sep = "\t", file = output_file)
            print("\t".join(line[:13]), "\t".join(cd), sep = "\t", file = output_b)
        else:
            se[ab] = ":".join(line)
    for line in ann_tidy:
        line = line.rstrip("\n").split("\t")
        if len(line) >1:
            ef = line[0] + ":" + line[1]
            if ef != ij:
                if ef in se.keys():
                    ghi = se[ef].split(":")
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
                    elif "stoploss" == consequence:
                        consequence = "stop_lost"
                    impact_ann = line[10]
                    gene_ann = line[11].split("-")
                    gene_ann = gene_ann[1]
                    coding_ann = line[13]
                    protein_ann = line[14].split(",")
                    protein_ann = protein_ann[0]
                    for i in range(len(line)):
                        if "0/1" in line[i] or "0/2" in line[i] or "0/3" in line[i]:
                            line[i] = "het"
                        elif "1/1" in line[i] or "2/2" in line[i] or "1/2" in line[i] or "1/3" in line[i] or "2/3" in line[i] or "2/1" in line[i]:
                            line[i] = "hom"
                        elif "./." in line[i]:
                            line[i] = "missing"
                        elif "0/0" in line[i]:
                            line[i] = "ref"
                    print("\t".join(ghi[:13]), SNP_ann, consequence, impact_ann,gene_ann, coding_ann, protein_ann, "\t".join(line[9:]), sep = "\t", file = output_file)
                    print("\t".join(ghi[:13]), SNP_ann, consequence, impact_ann,gene_ann, coding_ann, protein_ann, file = output_b, sep = "\t")
                    ij = ef

#find variants unique to breeds  
unique = []
u_b = []
for item in gb.keys():
    a = gb[item].split(":")
    if len(set(a)) <3:
        unique.append(item)
        u_b.append(a[1])

with open(path + "/lof_details.txt", "r") as input_file, open(path + "/unique_lof.txt", "w") as output_file, open(path + "/unique_lof_brief.txt", "w") as output_b:
    print("Breed\tgroup\tchrom\tpos\tref\talt\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein\tSNP_ann\tconsequence_ann\timpact_ann\tgene_ann\tcoding_ann\tprotein_ann", "\t".join(header[9:]), sep = "\t", file = output_file)
    print("Breed\tgroup\tchrom\tpos\tref\talt\tAC\tAF\tSNP\tconsequence\timpact\tgene\tcoding\tprotein\tSNP_ann\tconsequence_ann\timpact_ann\tgene_ann\tcoding_ann\tprotein_ann", sep = "\t", file = output_b)
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[1] + ":" + line[2]
        for i in range(len(unique)):
            if a == unique[i]:
                if u_b[i] != "Other":
                    print(u_b[i], "\t".join(line), file = output_file, sep ="\t")
                    print(u_b[i],"\t".join(line[:20]),file = output_b, sep = "\t")
