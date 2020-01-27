#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:59:25 2019

@author: durwa004
"""

directory = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/"

directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/"

count = 0
with open(directory + "/known_variants_locations.txt",encoding = "ISO-8859-1") as input_file, open(directory + "known_SNVs_tabix.sh", "w") as output_file:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        start = int(line[6]) - 50
        end = int(line[6]) + 50
        print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz ", line[5], ":", start, "-", end, " > ", line[1], "_", count, ".txt", sep = "", file = output_file)

#Do for QTLs (on MSI)    

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
        print("zgrep ", line[8], " GCA_000002305.1.refseq_chrs.vcf.gz > ", line[1], "_", line[0], ".txt", file = output_file, sep = "")
        rsid.append(line[8])

#Need to get a list of chrom:pos so that I can put into the remapping tool
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/known_qtls/split_files/"
with open("QTL_EC2_chrom_pos.txt", "w") as f:
    for filename in os.listdir(path + "/../"):
        if filename.endswith(".txt"):
            with open(path + "/" + filename) as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    print(line[0], ":", line[1], sep = "", file = f)
