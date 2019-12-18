#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:42:13 2019

@author: durwa004
"""

path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/"


#Goal = disease/horse/breed/genotype
#Get disease phenotypes
horse = []
breed = []
with open(path + "known_variants_for_analysis.txt", "r") as input_file, open(path + 
         "known_variants_with_breed.txt", "w") as f:
    print("Phenotype\thorse\tbreed\tgenotype\tdeleterious\tcausative", file = f)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "Phenotype" == line[0]:
            for i in range(len(line)):
                horse.append(line[i])
            next
        elif "NA" == line[0]:
            for i in range(len(line)):
                breed.append(line[i])
            next
        else:
            for i in range(len(line)):
                if line[i] == "1" or line[i] == "2":
                    print(line[1], horse[i], breed[i], line[i], line[2], line[3], file =f, sep = "\t")


        
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
            if "CC_" in a[1] or "Gait" in a[1] or "Curly" in a[1]:
                if "0/1" in line[i]:
                    print(1,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output_file)
                elif "1/1" in line[i]:
                    print(2,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output_file)
            elif a[1] == "ED":
                if "0/1" in line[i]:
                    print(1,"HERDA","\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output2)
                elif "1/1" in line[i]:
                    print(2,"HERDA","\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output2)
            elif a[1] == "ED_4":
                if "0/1" in line[i]:
                    print(1,"WFFS","\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output2)
                elif "1/1" in line[i]:
                    print(2,"WFFS","\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output2)
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

#Get number of known disease variants for each indidivual horse
with open(path + "known_variants_AFs.txt", "r"
          ) as input_file, open(path + "known_dz_genotypes_.txt", "w") as output2:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = phenotype[line[1]]
            a = a.split(":")
            if "CC_" in a[1] or "Gait" in a[1] or "Curly" in a[1]:
                next
            else:
                line1 = []
                for i in range(len(line)):
                    if "0/0" in line[i]:
                        line1.append("0")
                    elif "0/1" in line[i]:
                        line1.append("1")
                    elif "1/1" in line[i]:
                        line1.append("2")
        print("\t".join(line1), file = output2)
        
        
with open(path + "known_variants_AFs.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = phenotype[line[1]]
            a = a.split(":")
            if "IMM" in a[1]:
                if "0/1" in line[i]:
                    print("het:", vcf[i], ":", horse[vcf[i]])
                elif "1/1" in line[i]:
                    print("hom:", vcf[i], ":", horse[vcf[i]])
