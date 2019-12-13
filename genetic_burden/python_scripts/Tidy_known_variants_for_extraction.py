#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:59:25 2019

@author: durwa004
"""

directory = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/"

count = 0
with open(directory + "OMIA_known_SNVs.txt",encoding = "ISO-8859-1") as input_file, open(directory + 
         "known_SNVs_tabix.sh", "w") as output_file, open(directory + "known_SNVs.txt", "w") as output:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        start = int(line[6]) - 100
        end = int(line[7]) + 100
        print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz ", 
              line[2], ":", start, "-", end, " > ", line[0], "_", line[1], "_", count, ".txt", sep = "", file = output_file)
        print("\t".join(line), file = output)        


#Do for QTLs (on MSI)    
#Need to get disease for each QTL as well as the tabix code
import os
directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/known_qtls"
chrom_pos = {}
unmappable = {}
count = 0
for filename in os.listdir(directory):
    if filename.endswith("_report"):
        with open(directory + "/" + filename) as input_file:
            input_file.readline()
            for line in input_file:
                line = line.rstrip("\n").split("\t")
                if len(line) > 13:
                    if "First pass" in line[16]:
                        a = line[3] + ":" + line[8]
                        b = line[4] + ":" + line[13]
                        if a in chrom_pos.keys():
                            c = chrom_pos[a] + "," + b
                            chrom_pos[a] = c 
                        else:
                            chrom_pos[a] = b
                else:
                    print(line)
                    #a = line[3] + ":" + line[
        

with open(directory + "OMIA_known_SNVs.txt",encoding = "ISO-8859-1") as input_file, open(directory +
         "known_SNVs_tabix.sh", "w") as output_file, open(directory + "known_SNVs.txt", "w") as output:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        start = int(line[6]) - 100
        end = int(line[7]) + 100
        print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz ",
              line[2], ":", start, "-", end, " > ", line[0], "_", line[1], "_", count, ".txt", sep = "", file = output_file)
        print("\t".join(line), file = output)


