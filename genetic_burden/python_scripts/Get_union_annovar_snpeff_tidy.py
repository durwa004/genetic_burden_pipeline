import gzip
import os

#Goal: get union/intersect just based on whether annovar/snpeff called them coding
##Also get union/intersect based on whether annovar/snpeff called them high/moderate exactly
##Also get union/intersect based on whether annovar/snpeff called them high/moderate/low exactly

snpeff_dict = {}
with open("SnpEff_coding_tidy.txt","r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NW" in line[0]:
            next
        else:
            a = line[0] + ":" + line[1]
            snpeff_dict[a] = line[10]

len(set(snpeff_dict))
#Number high/mod/low
snpeff_impact = list(snpeff_dict.values())
print(snpeff_impact.count("HIGH"))
print(snpeff_impact.count("MODERATE"))
print(snpeff_impact.count("LOW"))

annovar_dict = {}
se_ann_exact = {}
with open("annovar_coding_tidy.txt","r") as input_file, open("snpeff_annovar_exact_intersect.txt","w") as output_file:
    input_file.readline()
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NW" in line[0]:
            next
        else:
            a = line[0] + ":" + line[1]
            annovar_dict[a] = line[10]
#Try and build in getting intersect to this as well - get exact intersect (if match for the impact)
            if a in snpeff_dict.keys():
                if snpeff_dict[a] == line[10]:
                    print(line[0], line[1], line[10], file = output_file)
                    se_ann_exact[a] = line[10]
se_ann_impact = list(se_ann_exact.values())
print(len(se_ann_impact))
print(se_ann_impact.count("HIGH"))
print(se_ann_impact.count("MODERATE"))
print(se_ann_impact.count("LOW"))

annovar_impact = list(annovar_dict.values())
print(len(annovar_impact))
print(annovar_impact.count("HIGH"))
print(annovar_impact.count("MODERATE"))
print(annovar_impact.count("LOW"))

#Extract just chrom/pos of the high/moderate impact variants
with open("snpeff_annovar_exact_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for key in se_ann_exact.keys():
        if se_ann_exact[key] == "HIGH" or se_ann_exact[key] == "MODERATE":
            a = key.split(":")
            print(a[0], a[1], se_ann_exact[key], sep = "\t", file = output_file)

#Get number of variants in each category for intersect of combined high/mod (not necessarily an exact match)
se_high_mod = {}
for key in snpeff_dict.keys():
    if snpeff_dict[key] == "HIGH" or snpeff_dict[key] == "MODERATE": 
        se_high_mod[key] = snpeff_dict[key]
print(len(se_high_mod))

ann_high_mod = {}
for key in annovar_dict.keys():
    if annovar_dict[key] == "HIGH" or annovar_dict[key] == "MODERATE":
        ann_high_mod[key] = annovar_dict[key]
print(len(ann_high_mod))

union = se_high_mod
union.update(ann_high_mod)
print(len(union))
union_impact = list(union.values())
print(union_impact.count("HIGH"))
print(union_impact.count("MODERATE"))

with open("snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for key in union.keys():
        a = key.split(":")
        print(a[0], a[1], union[key], sep = "\t", file = output_file)
            
#Get union and intersect of snpeff and annovar on position alone
intersect_pos = {}
for key in snpeff_dict.keys():
    if key in annovar_dict.keys():
        intersect_pos[key] = snpeff_dict[key]
print(len(intersect_pos))

union_dict = snpeff_dict
union_dict.update(annovar_dict)
print(len(union_dict))

