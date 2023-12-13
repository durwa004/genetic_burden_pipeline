#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:59:25 2019

@author: durwa004
"""

count = 0

with open("known_variants.txt",encoding = "ISO-8859-1") as input_file, open("known_variants_tabix.sh", "w") as output_file:
    input_file.readline()
    for line in input_file:
        count +=1
        line = line.rstrip("\n").split("\t")
        location = line[4] + ":" line[5] + "-" + line[6]
        print("tabix /home/durwa004/shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz ", location, " > ", line[0], ".txt", sep = "", file = output_file)
