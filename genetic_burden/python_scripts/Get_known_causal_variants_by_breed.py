#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:42:13 2019

@author: durwa004
"""

path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/"

#Get breed info
horse = {}
with open(path + "horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse[line[0]] = line[1]
horse['TWILIGHT'] = "TB"

#Get order of horse ids in vcf
vcf = []
with open(path + "horse_ids_vcf.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n")
        if line == "QUAL":
            vcf.append("AC")
        elif line == "FILTER":
            vcf.append("AF")
        elif line == "INFO" or line == "FORMAT":
            next
        else:
            vcf.append(line)

#Get disease phenotypes
phenotype = {}
with open(path + "causal_variants_with_chrom_pos.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[2] + ":" + line[1]
        phenotype[line[3]] = a
        
with open(path + "known_variants_AFs.txt", "r"
          ) as input_file, open(path + "known_CC_with_breeds.txt", "w"
          ) as output_file, open(path + "known_dz_with_breeds.txt", "w") as output2:
    print("Hom_het\tdisease\tCHROM\tPOS\tID\tREF\tALT\tAC\tAF\tbreed", file = output_file)
    print("Hom_het\tdisease\tCHROM\tPOS\tID\tREF\tALT\tAC\tAF\tbreed", file = output2)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = phenotype[line[1]]
            a = a.split(":")
            if "CC_" in a[1] or "Gait" in a[1]:
                if "0/1" in line[i]:
                    print(1,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output_file)
                elif "1/1" in line[i]:
                    print(2,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output_file)
            else:
                if "0/1" in line[i]:
                    print(1,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output2)
                elif "1/1" in line[i]:
                    print(2,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output2)

het_C = 0
hom_C = 0  
het_D = 0
hom_D = 0
with open(path + "known_variants_AFs.txt", "r"
          ) as input_file, open(path + "known_CC_with_breeds_hom_het.txt", "w"
          ) as output_file, open(path + "known_dz_with_breeds_hom_het.txt", "w") as output2:
    print("disease\tcount_het\tcount_hom", file = output_file)
    print("disease\tcount_het\tcount_hom", file = output2)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = phenotype[line[1]]
            a = a.split(":")
            if "CC_" in a[1] or "Gait" in a[1]:
                if "0/1" in line[i]:
                    het_C +=1
                elif "1/1" in line[i]:
                    hom_C +=1
            else:
                if "0/1" in line[i]:
                    het_D+=1
                elif "1/1" in line[i]:
                    hom_D+=1
        print(a[1],het_D,hom_D, sep = "\t", file = output2)
        print(a[1],het_C,hom_C, sep = "\t", file = output_file)
