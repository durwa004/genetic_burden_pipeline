import gzip
import os

#Goal: get union/intersect just based on whether annovar/snpeff called them coding
##Also get union/intersect based on whether annovar/snpeff called them high/moderate exactly
##Also get union/intersect based on whether annovar/snpeff called them high/moderate/low exactly

snpeff_chrom = []
snpeff_pos = []
snpeff_impact = []
snpeff_chrom_pos = []
snpeff_dict = {}
with open("SnpEff_coding_tidy.txt","r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        snpeff_chrom.append(line[0])
        snpeff_pos.append(line[1])
        snpeff_impact.append(line[10])
        a = line[0] + ":" + line[1]
        snpeff_chrom_pos.append(a)
        snpeff_dict[a] = line[10]

len(set(snpeff_chrom_pos))
#Number high/mod/low
snpeff_impact.count("HIGH")
snpeff_impact.count("MODERATE")
snpeff_impact.count("LOW")

annovar_chrom = []
annovar_pos = []
annovar_impact = []
annovar_chrom_pos = []
snpeff_annovar_chrom_pos = []
snpeff_annovar_impact = []

with open("annovar_coding_tidy.txt","r") as input_file, open("snpeff_annovar_exact_intersect.txt","w") as output_file:
    input_file.readline()
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        annovar_chrom.append(line[0])
        annovar_pos.append(line[1])
        annovar_impact.append(line[10])
        a = line[0] + ":" + line[1]
        annovar_chrom_pos.append(a)
#Try and build in getting intersect to this as well - get exact intersect (if match for the impact)
#Takes ~ 36 hours to run!
        for i,v in enumerate(snpeff_chrom_pos):
            if v == a:
                if snpeff_impact[i] == line[10]:
                    print(line[0], line[1], line[10], file = output_file)
                    snpeff_annovar_chrom_pos.append(a)
                    snpeff_annovar_impact.append(line[10])

#Extract just chrom/pos so that I can pull the intersect out of the snpeff file
with open("snpeff_annovar_exact_intersect_chrom_pos.txt", "w") as output_file:
    for i in range(len(snpeff_annovar_chrom_pos)):
        a = snpeff_annovar_chrom_pos[i].split(":")
        print(a[0],a[1], sep = "\t", file = output_file)

#Extract just chrom/pos of the high/moderate impact variants
with open("snpeff_annovar_exact_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for i in range(len(snpeff_annovar_chrom_pos)):
        if snpeff_annovar_impact[i] == "HIGH" or snpeff_annovar_impact[i] == "MODERATE":
            a = snpeff_annovar_chrom_pos[i].split(":")
            print(a[0],a[1], sep = "\t", file = output_file)
        else:
            continue

len(set(annovar_chrom_pos))
snpeff_annovar_chrom_pos = annovar_chrom_pos + snpeff_chrom_pos
len(set(snpeff_annovar_chrom_pos))

#Get union/intersect of snpeff and annovar coding variants:
union = list(set(snpeff_chrom_pos + annovar_chrom_pos))
intersect = list(set(snpeff_chrom_pos) & set(annovar_chrom_pos))
len(intersect)
#Get union/intersect if the variant is called high/moderate by both programs (do not have to agree on high/moderate
snpeff_high_mod_chrom_pos = []
snpeff_high_mod_impact = []
for i,v in enumerate(snpeff_chrom_pos):
    if snpeff_impact[i] == "HIGH" or snpeff_impact[i] == "MODERATE":
       snpeff_high_mod_chrom_pos.append(v)
       snpeff_high_mod_impact.append(snpeff_impact[i])
    else:
        next

annovar_high_mod_chrom_pos = []
annovar_high_mod_impact = []
for i,v in enumerate(annovar_chrom_pos):
    if annovar_impact[i] == "HIGH" or annovar_impact[i] == "MODERATE":
       annovar_high_mod_chrom_pos.append(v)
       annovar_high_mod_impact.append(annovar_impact[i])
    else:
       next

union_high_mod = list(set(snpeff_high_mod_chrom_pos + annovar_high_mod_chrom_pos))
intersect_high_mod = list(set(snpeff_high_mod_chrom_pos) & set(annovar_high_mod_chrom_pos))

len(set(union_high_mod))
len(set(intersect_high_mod))
#Get number of variants in each category for intersect of combined high/mod (not necessarily an exact match)
intersect_high_mod_impact = []
for i in range(len(intersect_high_mod)):
    for item in range(len(snpeff_high_mod_chrom_pos)):
        if intersect_high_mod[i] == snpeff_high_mod_chrom_pos[item]:
            intersect_high_mod_impact.append(snpeff_high_mod_impact[item])

intersect_high_mod_impact.count("HIGH")
intersect_high_mod_impact.count("MODERATE")

with open("snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for i in range(len(union_high_mod)):
        a = union_high_mod[i].split(":")
        print(a[0],a[1], sep = "\t", file = output_file)
