import gzip

snpeff = {}
with open("SnpEff.hml.txt", "r") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        if line[7] == "MODIFIER":
            next
        else:
            snpeff[a] = line[7]
se = list(snpeff.values())
print(se.count("HIGH"))
print(se.count("MODERATE"))
print(se.count("LOW"))
print(se.count("MODIFIER"))

VEP = {}
with open("VEP.hml.txt", "r") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        if line[7] == "MODIFIER":
            next
        else:
            VEP[a] = line[7]

ve = list(VEP.values())
print(ve.count("HIGH"))
print(ve.count("MODERATE"))
print(ve.count("LOW"))
print(ve.count("MODIFIER"))

intersect = {}
snpeff_only = {}
VEP_only = {}
for key in snpeff.keys():
    if key in VEP.keys():
        a = snpeff[key] + "_" + VEP[key]
        intersect[key] = a
    else:
        snpeff_only[key] = snpeff[key]

for key in VEP.keys():
    if key in snpeff.keys():
        if key in intersect.keys():
            next
        else:
            a = snpeff[key] + "_" + VEP[key] 
            intersect[key] = a
    else:
        VEP_only[key] = VEP[key]

se = list(snpeff_only.values())
print(se.count("HIGH"))
print(se.count("MODERATE"))
print(se.count("LOW"))

ve = list(VEP_only.values())
print(ve.count("HIGH"))
print(ve.count("MODERATE"))
print(ve.count("LOW"))

with open("SnpEff_VEP_intersect.txt", "w") as output_f:
    for key in intersect.keys():
        print(key, intersect[key], file = output_f)

##Transfer to laptop to create table of SnpEff/VEP overlap
##Get only H:H, M:H, and H:M variants
snpeff = {}

with open("SnpEff.hml.txt", "r") as input_f, open("SnpEff.intersect.txt", "w") as outputf:
#with open("SnpEff.hml.txt", "r") as input_f, open("SnpEff.mod.intersect.txt", "w") as outputf:
    a = input_f.readline()
    a = a.rstrip("\n").split("\t")
    alpha = ":".join(a[0:6])
    beta = ["SE." + x for x in a[6:]]
    snpeff[alpha] = ":".join(beta)
    print("Intersect\tKey", "\t".join(a), sep = "\t", file = outputf)
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        b = line[0] + ":" + line[1]
        for i,v in enumerate(line):
            if v == "":
                line[i] = "NA"
        if b in intersect.keys():
            if intersect[b] == "HIGH_HIGH" or intersect[b] == "MODERATE_HIGH" or intersect[b] == "HIGH_MODERATE":
#            if intersect[b] == "MODERATE_MODERATE" or intersect[b] == "MODERATE_LOW" or intersect[b] == "LOW_MODERATE":
                print(intersect[b], key, "\t".join(line), sep = "\t", file = outputf)
                c = ":".join(line[6:])
                snpeff[b] = c
            else:
                next

with open("VEP.hml.txt", "r") as input_f, open("VEP.intersect.txt", "w") as outputf, open("SnpEff.VEP.intersect.txt", "w") as outputf2:
#with open("VEP.hml.txt", "r") as input_f, open("VEP.mod.intersect.txt", "w") as outputf, open("SnpEff.VEP.mod.intersect.txt", "w") as outputf2: 
    for key in snpeff.keys():
        if "chrom" in key:
            aa = key.split(":")
            bb = snpeff[key].split(":")
    a = input_f.readline()
    a = a.rstrip("\n").split("\t")
    cc = ["VEP." + x for x in a[6:]]
    for i,v in enumerate(cc):
        bb.append(v)
    print("Intersect", "\t".join(aa), "\t".join(bb), sep = "\t", file = outputf2)
    print("Intersect\tKey", "\t".join(a), sep = "\t", file = outputf)
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        b = line[0] + ":" + line[1]
        for i,v in enumerate(line):
            if v == "":
                line[i] = "NA"
        if b in intersect.keys():
            if intersect[b] == "HIGH_HIGH" or intersect[b] == "MODERATE_HIGH" or intersect[b] == "HIGH_MODERATE":
            if intersect[b] == "MODERATE_MODERATE" or intersect[b] == "MODERATE_LOW" or intersect[b] == "LOW_MODERATE":
                print(intersect[b], key, "\t".join(line), sep = "\t", file = outputf)
                dd = snpeff[b].split(":")
                print(intersect[b], "\t".join(line[0:6]), "\t".join(dd), "\t".join(line[6:]), sep = "\t", file = outputf2)
   


#####################May not use

#Number high/mod/low
snpeff_impact = list(intersect.values())
print(snpeff_impact.count("HIGH_HIGH"))
print(snpeff_impact.count("HIGH"))
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
                lof_snpeff[a] = "LOF"
lof_se = list(lof_snpeff.values())
print(len(lof_se))

#Update the snpeff_dict values to LOF instead of high
for key in lof_snpeff.keys():
    snpeff_dict[key] = lof_snpeff[key]

#same for annvoar
annovar_dict = {}
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
                    lof_annovar[a] = "LOF"
lof_ann = list(lof_annovar.values())
print(len(lof_ann))

#Update the annovar_dict values to LOF instead of high
for key in lof_annovar.keys():
    annovar_dict[key] = lof_annovar[key]

#Get number of variants in each category for intersect of combined high/mod (not necessarily an exact match)
unique = 0
for key in snpeff_dict.keys():
    if key in annovar_dict.keys():    
        if snpeff_dict[key] == annovar_dict[key]:
            a = snpeff_dict[key] + ":" + annovar_dict[key]
            se_ann_union[key] = a
        else:
            a = snpeff_dict[key] + ":" + annovar_dict[key]
    else:
        a = snpeff_dict[key] + ":snpeff_only"
        se_ann_union[key] = a
        unique +=1

for key in annovar_dict.keys():
    if key in se_ann_union.keys():
        continue
    else:
        a = annovar_dict[key] + ":annovar_only"
        se_ann_union[key] = a
        unique +=1
print(len(se_ann_union))

se_ann_un= list(se_ann_union.values())
for i in list(set(se_ann_union.values())):
    print(i, ":", se_ann_un.count(i))

#Overall length of intersect = 
print(len(se_ann_un) - unique
#Get number of LOF variants in each category for intersect of combined high/mod (not necessarily an exact match)
lof_intersect = {}
for item in lof_snpeff.keys():
    if item in lof_annovar.keys():
         lof_intersect[item] = lof_snpeff[item]

print(len(lof_intersect))

#Get those that are lof/lof or lof/high or high/lof
ann_d = {}
for item in lof_annovar.keys():
    if item in lof_intersect.keys():
        next
    else:
        if item in snpeff_dict.keys():
            ann_d[item] = snpeff_dict[item]
        else:
            ann_d[item] = "missing"

se_d = {}
for item in lof_snpeff.keys():
    if item in lof_intersect.keys():
        next
    else:
        if item in annovar_dict.keys():
            se_d[item] = annovar_dict[item]
        else:
            se_d[item] = "missing"

#Only high impact in lof variants
lof_ann_only = list(ann_d.values())
for i in list(set(lof_ann_only)):
    print(i, ":", lof_ann_only.count(i))

lof_se_only = list(se_d.values())
for i in list(set(lof_se_only)):
    print(i, ":", lof_se_only.count(i))

#Extract just chrom/pos of the high/moderate impact variants
count =0
with open(path + "/gb_analysis/snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "w") as output_file:
    for key in se_ann_union.keys():
        if se_ann_union[key] == "HIGH:HIGH" or se_ann_union[key] == "HIGH:MODERATE" or se_ann_union[key] == "MODERATE:HIGH" or se_ann_union[key] == "LOF:LOF" or se_ann_union[key] == "LOF:HIGH" or se_ann_union[key] == "HIGH:LOF" or se_ann_union[key] == "LOF:MODERATE" or se_ann_union[key] == "MODERATE:LOF":
            a = key.split(":")
            b = se_ann_union[key].split(":")
            print(a[0], a[1], b[0], b[1], sep = "\t", file = output_file)
            count +=1

#Extract all lof variants called lof by both or lof by one and high by the other
count = 0
with open("lof/lof_combined_intersect_lof_high_chrom_pos.txt", "w") as output_file:
    for key in se_ann_union.keys():
        if se_ann_union[key] == "LOF:LOF" or se_ann_union[key] == "LOF:HIGH" or se_ann_union[key] == "HIGH:LOF":
            a = key.split(":")
            print(a[0], a[1], se_ann_union[key], sep = "\t", file = output_file)
            count +=1

