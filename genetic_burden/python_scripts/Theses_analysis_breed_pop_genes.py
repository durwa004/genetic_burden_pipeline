#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 05:01:00 2019

@author: durwa004
"""

import os

directory = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/"
#Common breed rare pop genes
gene = []
with open(directory + "common_breed_rare_pop_snpeff_info.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        gene.append(line[6])

len(gene)
unique = list(set(gene))
len(unique) #1478
with open(directory + "common_breed_rare_pop_unique_genes.txt","w") as output_file:
    for i in range(len(unique)):
        print(unique[i], file = output_file)

        
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol
original = {}
with open(directory + "common_breed_rare_pop_unique_genes_with_IDs.txt", 
          "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] == "-":
            next
        else:
            original[line[0]] = line[1]
len(original) #1,218

with open(directory + "common_breed_rare_pop_snpeff_info_with_genes.txt", 
             "w") as output_file, open(directory + "common_breed_rare_pop_snpeff_info.txt", 
                "r") as input2:
    input2.readline()
    for line in input2:
        line = line.rstrip("\n").split("\t")
        b = line[6].split(".")
        if b[0] in original.keys():
            print("\t".join(line[0:6]),original[b[0]],"\t".join(line[7:]),sep="\t",file = output_file)
        else:
            print("\t".join(line), file = output_file)
        
with open(directory + "common_breed_rare_pop_snpeff_info_with_genes.txt", 
          "r") as input_file, open(directory + "common_breed_rare_pop_genes_coding.txt", 
             "w") as coding, open(directory + "common_breed_rare_pop_genes_modifier.txt", 
                "w") as modifier:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[5] == "MODIFIER":
            print(line[6], file = modifier)
        else:
            print(line[6], file = coding)
#Check enrichment: https://amp.pharm.mssm.edu/Enrichr/enrich
            
            
directory = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/"
#Rare breed common pop genes
gene = []
with open(directory + "rare_breed_common_pop_snpeff_info.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        gene.append(line[6])

len(gene)
unique = list(set(gene))
len(unique) #26,007

with open(directory + "rare_breed_common_pop_unique_genes1.txt",
          "w") as output_file, open(directory + "rare_breed_common_pop_unique_genes2.txt", 
             "w") as output2, open(directory + "rare_breed_common_pop_unique_genes3.txt",
             "w") as output3, open(directory + "rare_breed_common_pop_unique_genes4.txt", 
             "w") as output4, open(directory + "rare_breed_common_pop_unique_genes5.txt", "w") as output5:
    print("\n".join(unique[0:5000]), file = output_file)
    print("\n".join(unique[5000:10000]), file = output2)
    print("\n".join(unique[10000:15000]), file = output3)
    print("\n".join(unique[15000:20000]), file = output4)
    print("\n".join(unique[20000:]), file = output5)
        
#Had to split this list into 5 to try and avoid missing entries
        
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol
        
original = {}
with open(directory + "rare_breed_common_pop_unique_genes_with_IDs.txt", 
          "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] == "-":
            next
        else:
            original[line[0]] = line[1]
len(original) #19,348

with open(directory + "rare_breed_common_pop_snpeff_info_with_genes.txt", 
             "w") as output_file, open(directory + "rare_breed_common_pop_snpeff_info.txt", 
                "r") as input2:
    input2.readline()
    for line in input2:
        line = line.rstrip("\n").split("\t")
        b = line[6].split(".")
        if b[0] in original.keys():
            print("\t".join(line[0:6]),original[b[0]],"\t".join(line[7:]),sep="\t",file = output_file)
        else:
            print("\t".join(line),file=output_file)
        
with open(directory + "rare_breed_common_pop_snpeff_info_with_genes.txt", 
          "r") as input_file, open(directory + "rare_breed_common_pop_genes_coding.txt", 
             "w") as coding, open(directory + "rare_breed_common_pop_genes_modifier.txt", 
                "w") as modifier:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[5] == "MODIFIER":
            print(line[6], file = modifier)
        else:
            print(line[6], file = coding)
#Check enrichment: https://amp.pharm.mssm.edu/Enrichr/enrich     

#Which genes have high impact variants in them? 
gene_high = []
gene_moderate = []
gene_low = []
gene_modifier = []

with open(directory + "rare_breed_common_pop_snpeff_info_with_genes.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[5] == "HIGH":
            gene_high.append(line[6])
        elif line[5] == "MODERATE":
            gene_moderate.append(line[6])
        elif line[5] == "LOW":
            gene_low.append(line[6])
        elif line[5] == "MODIFIER":
            gene_modifier.append(line[6])
        else:
            print(line[6])

