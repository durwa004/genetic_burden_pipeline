#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:59:25 2019

@author: durwa004
"""

directory = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/"

count = 0
with open(directory + "OMIA_known_SNVs_191211.txt",encoding = "ISO-8859-1") as input_file, open(directory + 
         "known_SNVs_tabix.sh", "w") as output_file, open(directory + "known_SNVs.txt", "w") as output:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        start = int(line[3]) - 100
        end = int(line[4]) + 100
        print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz ", 
              line[2], ":", start, "-", end, " > ", line[0], "_", line[1], "_", count, ".txt", sep = "", file = output_file)
        print("\t".join(line), file = output)        

#Do for QTLs (on MSI)    
#Need to get disease for each QTL as well as the tabix code
import os
import numpy as np
directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/known_qtls/split_files/"
chrom_pos = {}
chrom_pos_2 = {}
unmappable = {}
count_m = 0
count_u = 0
count_2 = 0
count = 0
for filename in os.listdir(directory):
    if filename.endswith("_report"):
        with open(directory + "/" + filename) as input_file:
            input_file.readline()
            for line in input_file:
                line = line.rstrip("\n").split("\t")
                count +=1
                if len(line) > 13:
                    if "First Pass" in line[16]:
                        count_m +=1
                        a = line[3] + ":" + line[8]
                        b = line[4] + ":" + line[13]
                        if a in chrom_pos.keys():
                            c = chrom_pos[a] + "," + b
                            chrom_pos[a] = c 
                        else:
                            chrom_pos[a] = b
                    elif "Second Pass" in line[16]:
                        count_2 +=1
                        a = line[3] + ":" + line[8]
                        b = line[4] + ":" + line[13]
                        if a in chrom_pos_2.keys():
                            c = chrom_pos_2[a] + "," + b
                            chrom_pos_2[a] = c
                        else:
                            chrom_pos_2[a] = b
                else:
                    unmappable[line[3]] = line[5]
                    count_u +=1        

np.setdiff1d(list(chrom_pos_2.keys()), list(chrom_pos.keys()))
chrom_pos["NC_009157:1368081"] = "NC_009172.3:1876173"
chrom_pos["NC_009157:1388861"] = "NC_009172.3:1876173"
chrom_pos["NC_009157:1427118"] = 'NC_009172.3:1876173'

#Need to figure out what dz they were associated with before I run this


count = 0
with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/QTLs_remapped.sh", "w") as f:
    for key, value in chrom_pos.items():
        a = value.split(",")
        if len(a) >1:
            count +=1
        else:
            print(

with open(directory + "OMIA_known_SNVs.txt",encoding = "ISO-8859-1") as input_file, open(directory +
         "known_SNVs_tabix.sh", "w") as output_file, open(directory + "known_SNVs.txt", "w") as output:

#Get QTLs to pull out using grep
count = 0
rsid = []
with open(directory + "animalgenomeQTL_for_extraction.txt",encoding = "ISO-8859-1") as input_file, open(directory + 
         "grep_get_rs_chrom_pos.sh", "w") as output_file:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        start = int(line[6]) - 100
        end = int(line[7]) + 100
        print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz ",
              line[2], ":", start, "-", end, " > ", line[0], "_", line[1], "_", count, ".txt", sep = "", file = output_file)
        print("\t".join(line), file = output)


=======
        print("zgrep ", line[8], " GCA_000002305.1.refseq_chrs.vcf.gz > ", line[1], "_", line[0], ".txt", file = output_file, sep = "")
        rsid.append(line[8])
<<<<<<< HEAD

#Need to get a list of chrom:pos so that I can put into the remapping tool
        
=======
 
#Need to get a list of chrom:pos so that I can put into the remapping tool
>>>>>>> 882e64fd6ea0b82f9567aebbb3714ab63edbef84
>>>>>>> 239cb49afcf71a214d7b01526853ce523d3c2662
