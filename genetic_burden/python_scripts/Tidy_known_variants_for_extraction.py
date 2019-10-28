#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:59:25 2019

@author: durwa004
"""

directory = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/"

count = 0
with open(directory + "known_causal_SNVs.txt",encoding = "ISO-8859-1") as input_file, open(directory + 
         "known_SNVs_tabix.sh", "w") as output_file, open(directory + "known_SNVs.txt", "w") as output:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        start = int(line[3]) - 100
        end = int(line[4]) + 100
        print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/thesis_intersect.vcf.gz ", 
              line[2], ":", start, "-", end, " > ", line[0], "_", line[1], "_", count, ".txt", sep = "", file = output_file)
        print("\t".join(line), file = output)        
    