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

#Get dictionary of just lof alleles
lof_snpeff = {} 
with open(path + "SnpEff/SnpEff_coding_tidy.txt","r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NW" in line[0]:
            next
        else:
            if line[15] == "y":
                a = line[0] + ":" + line[1]
                lof_snpeff[a] = line[10]
lof_se = list(lof_snpeff.values())
print(len(lof_se))
print(lof_se.count("HIGH"))
print(lof_se.count("MODERATE"))
print(lof_se.count("LOW"))


#same for annvoar
annovar_dict = {}
se_ann_exact = {}
with open(path + "annovar/annovar_coding_tidy.txt","r") as input_file:
    input_file.readline()
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NW" in line[0]:
            pass 
        else:
            a = line[0] + ":" + line[1]
            annovar_dict[a] = line[10]

annovar_impact = list(annovar_dict.values())
print(len(annovar_impact))
print(annovar_impact.count("HIGH"))
print(annovar_impact.count("MODERATE"))
print(annovar_impact.count("LOW"))

lof_annovar = {}
with open(path + "annovar/annovar_coding_tidy_lof.txt","r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if len(line) >1:
            if "NW" in line[0]:
                pass
            else:
                if line[15] == "y":
                    a = line[0] + ":" + line[1]
                    lof_annovar[a] = line[10]
lof_ann = list(lof_annovar.values())
print(len(lof_ann))
print(lof_ann.count("HIGH"))
print(lof_ann.count("MODERATE"))
print(lof_ann.count("LOW"))


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

#Get number of LOF variants in each category for intersect of combined high/mod (not necessarily an exact match)
lof_intersect = set(lof_annovar.keys()) & set(lof_snpeff.keys())
ann_only = set(lof_annovar.keys()) - set(lof_snpeff.keys())
snpeff_only = set(lof_snpeff.keys()) - set(lof_annovar.keys())

print(len(lof_intersect))
print(len(ann_only))
print(len(snpeff_only))

ann_only = list(ann_only)
snpeff_only = list(snpeff_only)
lof_intersect = list(lof_intersect)
#All lof variants have a high impact = intersect is all HIGH:HIGH
ann_d = {}
for item in range(len(ann_only)):
    if ann_only[item] in snpeff_dict.keys():
        ann_d[ann_only[item]] = snpeff_dict[ann_only[item]]
    else:
        ann_d[ann_only[item]] = "missing"

se_d = {}
for item in range(len(snpeff_only)):
    if snpeff_only[item] in annovar_dict.keys():
        se_d[snpeff_only[item]] = annovar_dict[snpeff_only[item]]
    else:
        se_d[snpeff_only[item]] = "missing"

#Only high impact in lof variants
lof_ann_only = list(ann_d.values())
print(lof_ann_only.count("HIGH"))
print(lof_ann_only.count("MODERATE"))
print(lof_ann_only.count("missing")) 

lof_se_only = list(se_d.values())
print(lof_se_only.count("HIGH"))
print(lof_se_only.count("MODERATE"))
print(lof_se_only.count("LOW"))
print(lof_se_only.count("UNKNOWN"))
print(lof_se_only.count("missing"))

lof_annovar.update(lof_snpeff)

#Extract just chrom/pos of the high/moderate impact variants
with open("snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for key in se_ann_intersect.keys():
        if se_ann_intersect[key] == "HIGH:HIGH" or se_ann_intersect[key] == "HIGH:MODERATE" or se_ann_intersect[key] == "MODERATE:HIGH":
            a = key.split(":")
            b = se_ann_intersect[key].split(":")
            print(a[0], a[1], b[0], b[1], sep = "\t", file = output_file)

#Extract all lof variants called lof by both or lof by one and high by the other
with open("../lof/lof_combined_intersect_lof_high_chrom_pos.txt", "w") as output_file:
    for i in range(len(lof_intersect)):
        a = lof_intersect[i].split(":")
        print(a[0], a[1], "intersect", sep = "\t", file = output_file)
    for key in ann_d.keys():
        if ann_d[key] == "HIGH":
            b = key.split(":")
            print(b[0], b[1], "ann_only", file = output_file, sep = "\t")
    for key in se_d.keys():
        if se_d[key] == "HIGH":
            b = key.split(":")
            print(b[0], b[1], "se_only", file = output_file, sep = "\t")


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
