#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:11:06 2019

@author: durwa004
"""

import os
import gzip

path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/no_homozygotes/"

#Get list of horse ids in order of vcf.
het = {}
hom = {}
missing = {}
header = []
with gzip.open(path + "../thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                het[line[i]] = 0
                hom[line[i]] = 0
                missing[line[i]] = 0
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

header1 = header[9:]

gb = {}
with gzip.open(path + "thesis_intersect_no_homozygotes.vcf.gz", "rt") as input_file, open(path + "/no_homozygotes_details.txt", "w") as output_file, open(path + "/lof_variants_no_homozygotes.txt", "w") as lof_file:
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = output_file)
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = lof_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t") 
        if "#" in line[0]:
            input_file.readline()
        else:
            c_p = line[0] + ":" + line[1]
            ab = line[7].split(";")
            cd = ab[1].split("AF=")
            AF = cd[1]
            if "," in AF:
                ef = AF.split(",")
                AF = ef[0]
                MA = "y"
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
            if "frameshift" in consequence or "start_lost" in consequence or "stop_gained" in consequence or "stop_lost" in consequence:
                lof = "y"
                print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = lof_file)
            else:
                lof = "n"
            print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = output_file)
            line1 = line[9:]
            for i in range(len(line1)):
                if "0/1" in line1[i] or "1/1" in line1[i] or "0/2" in line1[i] or "1/2" in line1[i] or "2/2" in line1[i] or "2/1" in line1[i] or "1/3" in line1[i] or "0/3" in line1[i] or "2/3" in line1[i]:
                     de = gb[c_p] + ":" + horse_breed[header1[i]]
                     gb[c_p] = de  

#find variants unique to breeds  
unique = []
u_b = []
for item in gb.keys():
    a = gb[item].split(":")
    if len(set(a)) <3:
        unique.append(item)
        u_b.append(a[1])

with open(path + "/no_homozygotes_details.txt", "r") as input_file, open(path + "/unique_no_homozygotes.txt", "w") as output_file:
    print("Breed\tCHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = output_file)
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        for i in range(len(unique)):
            if a == unique[i]:
                print(u_b[i], "\t".join(line), file = output_file, sep ="\t")


#All homozygotes
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/"

gb = {}
with gzip.open(path + "thesis_intersect_all_homozygotes.vcf.gz", "rt") as input_file, open(path + "/all_homozygotes_details.txt", "w") as output_file, open(path + "/lof_variants_all_homozygotes.txt", "w") as lof_file:
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = output_file)
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = lof_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            input_file.readline()
        else:
            c_p = line[0] + ":" + line[1]
            ab = line[7].split(";")
            cd = ab[1].split("AF=")
            AF = cd[1]
            if "," in AF:
                ef = AF.split(",")
                AF = ef[0]
                MA = "y"
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
            if "frameshift" in consequence or "start_lost" in consequence or "stop_gained" in consequence or "stop_lost" in consequence:
                lof = "y"
                print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = lof_file)
            else:
                lof = "n"
            print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = output_file)
            line1 = line[9:]
            for i in range(len(line1)):
                if "0/1" in line1[i] or "1/1" in line1[i] or "0/2" in line1[i] or "1/2" in line1[i] or "2/2" in line1[i] or "2/1" in line1[i] or "1/3" in line1[i] or "0/3" in line1[i] or "2/3" in line1[i]:
                     de = gb[c_p] + ":" + horse_breed[header1[i]]
                     gb[c_p] = de

#All homozygotes/heterozygotes
#Also want a separate file that excludes variants homozygous in all for above comparison
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/"

gb = {}
#with gzip.open(path + "thesis_intersect_variants_present_in_all.vcf.gz", "rt") as input_file, open(path + "/variants_present_in_all_details.txt", "w") as output_file, open(path + "/lof_variants_present_in_all.txt", "w") as lof_file, open(path + "/variants_present_in_all_not_homozygous_in_all.txt", "w") as non_hom_file:
#    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = output_file)
#    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = lof_file)
with gzip.open(path + "thesis_intersect_variants_present_in_all.vcf.gz", "rt") as input_file,  open(path + "/variants_present_in_all_not_homozygous_in_all.txt", "w") as non_hom_file:
    print("CHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = non_hom_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            input_file.readline()
        else:
            c_p = line[0] + ":" + line[1]
            ab = line[7].split(";")
            cd = ab[1].split("AF=")
            AF = cd[1]
            if "," in AF:
                ef = AF.split(",")
                AF = ef[0]
                MA = "y"
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
            if "frameshift" in consequence or "start_lost" in consequence or "stop_gained" in consequence or "stop_lost" in consequence:
                lof = "y"
#                print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = lof_file)
            else:
                lof = "n"
#            print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = output_file)
            line1 = line[9:]
            for i in range(len(line1)):
                if "0/1" in line1[i] or "1/1" in line1[i] or "0/2" in line1[i] or "1/2" in line1[i] or "2/2" in line1[i] or "2/1" in line1[i] or "1/3" in line1[i] or "0/3" in line1[i] or "2/3" in line1[i]:
                     de = gb[c_p] + ":" + horse_breed[header1[i]]
                     gb[c_p] = de
            count = 0
            for i in range(len(line1)):
                if "0/1" in line1[i] or "0/2" in line1[i] or "0/3" in line1[i]:
                    count +=1
            if count >0:
                print(line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = non_hom_file) 
