import gzip
import os

#Goal: get union/intersect just based on whether annovar/snpeff called them coding
##Also get union/intersect based on whether annovar/snpeff called them high/moderate exactly
##Also get union/intersect based on whether annovar/snpeff called them high/moderate/low exactly
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/"

snpeff_dict = {}
with open(path + "SnpEff/SnpEff_coding_tidy.txt","r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NW" in line[0]:
            next
        else:
            a = line[0] + ":" + line[1]
            snpeff_dict[a] = line[10]

print(len(set(snpeff_dict)))
#Number high/mod/low
snpeff_impact = list(snpeff_dict.values())
print(snpeff_impact.count("HIGH"))
print(snpeff_impact.count("MODERATE"))
print(snpeff_impact.count("LOW"))

annovar_dict = {}
se_ann_exact = {}
with open(path + "annovar/annovar_coding_tidy.txt","r") as input_file:
    input_file.readline()
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NW" in line[0]:
            next
        else:
            a = line[0] + ":" + line[1]
            annovar_dict[a] = line[10]

annovar_impact = list(annovar_dict.values())
print(len(annovar_impact))
print(annovar_impact.count("HIGH"))
print(annovar_impact.count("MODERATE"))
print(annovar_impact.count("LOW"))


#Get number of variants in each category for intersect of combined high/mod (not necessarily an exact match)
se_ann_union = {}
se_ann_intersect = {}

for key in snpeff_dict.keys():
    if key in annovar_dict.keys():
        a = snpeff_dict[key] + ":" + annovar_dict[key]
        se_ann_intersect[key] = a
        se_ann_union[key] = a
    else:
        a = snpeff_dict[key] + ":snpeff_only"
        se_ann_union[key] = a
print(len(se_ann_intersect))
print(len(se_ann_union))

for key in annovar_dict.keys():
    if key in snpeff_dict.keys():
        continue
    else:
        a = annovar_dict[key] + ":annovar_only"
        se_ann_union[key] = a
print(len(se_ann_union))

se_ann_int= list(se_ann_intersect.values())
print(se_ann_int.count("HIGH:HIGH"))
print(se_ann_int.count("HIGH:MODERATE"))
print(se_ann_int.count("MODERATE:HIGH"))
print(se_ann_int.count("MODERATE:MODERATE"))
print(se_ann_int.count("HIGH:LOW"))
print(se_ann_int.count("MODERATE:LOW"))
print(se_ann_int.count("LOW:LOW"))
print(se_ann_int.count("LOW:HIGH"))
print(se_ann_int.count("LOW:MODERATE"))

se_ann_un= list(se_ann_union.values())
print(se_ann_un.count("HIGH:HIGH"))
print(se_ann_un.count("HIGH:MODERATE"))
print(se_ann_un.count("MODERATE:HIGH"))
print(se_ann_un.count("MODERATE:MODERATE"))
print(se_ann_un.count("HIGH:LOW"))
print(se_ann_un.count("MODERATE:LOW"))
print(se_ann_un.count("LOW:LOW"))
print(se_ann_un.count("LOW:HIGH"))
print(se_ann_un.count("LOW:MODERATE"))
print(se_ann_un.count("HIGH:snpeff_only"))
print(se_ann_un.count("MODERATE:snpeff_only"))
print(se_ann_un.count("LOW:snpeff_only"))
print(se_ann_un.count("HIGH:annovar_only"))
print(se_ann_un.count("MODERATE:annovar_only"))
print(se_ann_un.count("LOW:annovar_only"))

#Extract just chrom/pos of the high/moderate impact variants
with open("snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for key in se_ann_intersect.keys():
        if se_ann_intersect[key] == "HIGH:HIGH" or se_ann_intersect[key] == "HIGH:MODERATE" or se_ann_intersect[key] == "MODERATE:HIGH":
            a = key.split(":")
            b = se_ann_intersect[key].split(":")
            print(a[0], a[1], b[0], b[1], sep = "\t", file = output_file)

#Intersect of variants
with open("snpeff_annovar_intersect.txt", "w") as output_file:
    for key in se_ann_intersect.keys():
        a = key.split(":")
        b = se_ann_intersect[key].split(":")
        print(a[0], a[1], b[0], b[1], sep = "\t", file = output_file)

#Union of variants
with open("snpeff_annovar_union.txt", "w") as output_file:
    for key in se_ann_union.keys():
        a = key.split(":")
        b = se_ann_union[key].split(":")
        print(a[0], a[1], b[0], b[1], sep = "\t", file = output_file)
