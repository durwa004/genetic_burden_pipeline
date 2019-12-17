#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:31:36 2019

@author: durwa004
"""

path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/"
genes = []
with open(path + "/breed_rare_other_breed_common.txt", "r") as input_file, open(path + 
         "/breed_rare_other_breed_common_genes.txt", "w") as f:
    for line in input_file:
        input_file.readline()
        line = line.rstrip("\n").split()
        genes.append(line[0])
    genes = list(set(genes))
    for i in range(len(genes)):
        print(genes[i], file = f)
#Use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ids
        #RefSeq mRNA accession to gene symbol

#Get breeds that most/least commonly have variant discrepancies
#Need to add the rare and common variants together for each breed difference
breed_disc = {}
breeds = []
with open(path + "/breed_differences.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        breed_disc[line[0]] = line[1]
        a = line[0].split("_rare_")
        breeds.append(a[0])
breeds = list(set(breeds))

br_dict = {}
for i in range(len(breeds)):
    for key,value in breed_disc.items():
        a = key.split("_rare_")
        b = a[1].split("_common")
        c = a[0] + ":" + b[0]
        e = b[0] + ":" + a[0]
        if c in br_dict.keys():
            d = int(br_dict[c]) + int(value)
            br_dict[c] = d
        elif e in br_dict.keys():
            d = int(br_dict[e]) + int(value)
            br_dict[e] = d
        else:
            br_dict[c] = int(value)

min_c = 100000000
max_c = 0
av_c = 0
count = 0
for key,value in br_dict.items():
    count +=1
    if int(value) < int(min_c):
        min_c = value
        lowest = key
    elif int(value) > int(max_c):
        max_c = value
        most = key
    av_c += int(value)

#Mean number of variant discrepancies between breeds
print(av_c/count)
#Smallest number of variant discrepancies between breeds
print(min_c)
print(lowest)
#Biggest number of variant discrepancies between breeds
print(max_c)
print(most)

#Reimport the gene names and create a dictionary
gene_ids = {}
with open("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/breed_rare_other_breed_common_genes_ids.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] == "-":
            next
        else:
            gene_ids[line[0]] = line[1]
len(gene_ids) #9,908

with open(path + "/breed_breed_disc_gene_list.txt" ,"w") as f:
    for key in gene_ids.keys():
        print(gene_ids[key], file = f)