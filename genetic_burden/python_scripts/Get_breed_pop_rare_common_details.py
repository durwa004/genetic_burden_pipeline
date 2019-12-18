#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 13:31:36 2019

@author: durwa004
"""

path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/"
#############################Breed common, rare pop
genes = []
with open(path + "/breed_common_pop.txt", "r") as input_file, open(path + 
         "/breed_common_pop_rare_genes.txt", "w") as f:
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
with open(path + "/breed_common_differences.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        breed_disc[line[0]] = line[1]

min_c = 100000000
max_c = 0
av_c = 0
count = 0
for key,value in breed_disc.items():
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
with open(path + "/breed_common_rare_pop_gene_ids.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] == "-":
            next
        else:
            gene_ids[line[0]] = line[1]
len(gene_ids) #9,908

with open(path + "/breed_common_pop_rare_gene_list.txt" ,"w") as f:
    for key in gene_ids.keys():
        print(gene_ids[key], file = f)
        
######################################Breed rare, common pop

path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/"
genes = []
with open(path + "/breed_common_pop.txt", "r") as input_file, open(path + 
         "/breed_common_pop_rare_genes.txt", "w") as f:
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
with open(path + "/breed_common_differences.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        breed_disc[line[0]] = line[1]

min_c = 100000000
max_c = 0
av_c = 0
count = 0
for key,value in breed_disc.items():
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
with open(path + "/breed_common_rare_pop_gene_ids.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] == "-":
            next
        else:
            gene_ids[line[0]] = line[1]
len(gene_ids) #9,908

with open(path + "/breed_rare_pop_common_gene_list.txt" ,"w") as f:
    for key in gene_ids.keys():
        print(gene_ids[key], file = f)
        
#####################################unique breed

path = "/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/"
genes = []
with open(path + "/breed_unique_pop.txt", "r") as input_file, open(path + 
         "/breed_unique_pop_genes.txt", "w") as f:
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
with open(path + "/breed_unique_differences.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        breed_disc[line[0]] = line[1]

min_c = 100000000
max_c = 0
av_c = 0
count = 0
for key,value in breed_disc.items():
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
with open(path + "/breed_unique_gene_ids.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] == "-":
            next
        else:
            gene_ids[line[0]] = line[1]
len(gene_ids) #9,908

with open(path + "/breed_unique_gene_list.txt" ,"w") as f:
    for key in gene_ids.keys():
        print(gene_ids[key], file = f)