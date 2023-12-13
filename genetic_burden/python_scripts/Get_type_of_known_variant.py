import os
import gzip


#Get details of frequency of known variants

#Get list of horse ids in order of vcf.
header = []
with gzip.open("../../shared/PopulationVCF/joint_genotype_combined.goldenPath.vep.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
            break

header[0] = "CHROM"
#Get variant details for each horse
count = len(header)
details = []
for i in range(count):
    details.append("NA")

with open("known_variants_genotypes.txt", "w") as output2:
    print("Phenotype\tAC\tAF\timpact\tconsequence", "\t".join(header), sep = "\t", file = output2)
    for filename in os.listdir("known_variants/"):
        if filename.endswith(".txt"):
            with open("known_variants/" + filename, "r") as input_file:
                if os.stat("known_variants/" + filename).st_size !=0:
                    for line in input_file:
                        line = line.rstrip("\n").split("\t")
                        filename1 = filename.split(".txt")
                        print(filename1[0], "\t".join(line),sep = "\t", file = output2)
                else:
                    filename1 = filename.split(".txt")
                    print(filename1[0], "\t".join(details), sep = "\t", file = output2)



with open("known_variants_genotypes.txt") as input_file, open("known_variants_tidy.txt", "w") as output_file:
    header = input_file.readline()
    header = header.rstrip("\n").split("\t")
    print("phenotype\tchrom\tpos\tref\talt\tAC\tAF\timpact\tconsequence\tgene_name\tgene_ID", "\t".join(header[10:]), sep = "\t", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        chr_pos = line[1] + ":" + line[2] + ":" + line[4] #chr:pos:ref
        horses = line[10:]
        for i,v in enumerate(horses):
            alpha = v.split(":")
            horses[i] = alpha[0]
        if "," in line[5]: 
            al = line[5].split(",")
            ab = line[8].split(";")
            bc = ab[0].split("AC=")
            bc = bc[1].split(",")
            cd = ab[1].split("AF=")
            cd = cd[1].split(",")
            for i in range(len(al)):
                xv = al[i] + ":" + bc[i] + ":" + cd[i]  #alt:AC:AF
                ef = line[8].split("CSQ=")
                ef = ef[1].split(";")
                ef = ef[0].split(",")
                if al[i] == "*":
                    xv1 = xv.split(":")
                    print(line[0], "\t".join(chr_pos1), "\t".join(xv1), "NA\tNA\tNA\tNA", "\t".join(horses), sep = "\t", file = output_file)
                for item1, value1 in enumerate(ef):
                    fg = value1.split("|")
                    if al[i] == fg[0]:
                        plant = []
                        for item, value in enumerate(fg[1:]):
                            if value == "":
                                plant.append("NA")
                            else:
                                plant.append(value)
                        chr_pos1 = chr_pos.split(":")
                        xv1 = xv.split(":")
                        print(line[0], "\t".join(chr_pos1), "\t".join(xv1), plant[1], plant[0], plant[2], plant[3], "\t".join(horses), sep = "\t", file = output_file)
                        break
                    else:
                        next
        else:
            al = line[5].split(",")
            ab = line[8].split(";")
            bc = ab[0].split("AC=")
            if ab[0] == "NA":
                bc = "NA"
                cd = "NA"
                plant = "NA:" * 10
                plant = plant.split(":")
            else:
                bc = bc[1]
                cd = ab[1].split("AF=")
                cd = cd[1]
                ef = line[8].split("CSQ=")
                ef = ef[1].split(";")
                ef = ef[0].split(",")
                fg = ef[0].split("|")
                plant = []
                for item, value in enumerate(fg[1:]):
                    if value == "":
                        plant.append("NA")
                    else:
                        plant.append(value)
            chr_pos1 = chr_pos.split(":")
            print(line[0],"\t".join(chr_pos1), line[5], bc, cd, plant[1], plant[0], plant[2], plant[3], "\t".join(horses), sep = "\t", file = output_file)

#Convert to table for R
horse = {}
with open("horse_genomes_breeds_tidy.txt", "r") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        if line[0] in horse.keys():
            next
        else:
            horse[line[0]] = line[1]

#Multiallelic sites are all first alternate except Lordosis (alt 2), mel3 (alt 2)
phenotype = {}
with open("known_variants_tidy.txt", "r") as input_f:
    header = input_f.readline()
    header = header.rstrip("\n").split("\t")
    header = header[11:]
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        if line[1] == "NA":
            phenotype[line[0]] = "NA"
        else:
            phenotype[line[0]] = {}
            breeds = list(set(horse.values()))
            for i, v in enumerate(breeds):
                phenotype[line[0]][v] = {}
                phenotype[line[0]][v]["het"] = 0
                phenotype[line[0]][v]["hom"] = 0
                phenotype[line[0]][v]["miss"] = 0
                phenotype[line[0]][v]["hom_WT"] = 0
            if line[0] == "Lord" or line[0] == "mel3":
                brief = line[11:]
                for i1, v1 in enumerate(brief):
                    if v1 == "0/0" or v1 == "0|0":
                        a = int(phenotype[line[0]][horse[header[i1]]]["hom_WT"]) + 1
                        phenotype[line[0]][horse[header[i1]]]["hom_WT"] = a
                    elif v1 == "0/2" or v1 == "0|2" or v1 == "1/2" or v1 == "1|2":
                        a = int(phenotype[line[0]][horse[header[i1]]]["het"]) + 1
                        phenotype[line[0]][horse[header[i1]]]["het"] = a
                    elif v1 == "2/2" or v1 == "2|2":
                        a = int(phenotype[line[0]][horse[header[i1]]]["hom"]) + 1
                        phenotype[line[0]][horse[header[i1]]]["hom"] = a
                    elif v1 == "./.":
                        a = int(phenotype[line[0]][horse[header[i1]]]["miss"]) + 1
                        phenotype[line[0]][horse[header[i1]]]["miss"] = a
            else:
                brief = line[11:]
                for i2, v2 in enumerate(brief):
                    if v2 == "0/0" or v2 == "0|0":
                        a = int(phenotype[line[0]][horse[header[i2]]]["hom_WT"]) + 1
                        phenotype[line[0]][horse[header[i2]]]["hom_WT"] = a
                    elif v2 == "0/1" or v2 == "0|1":
                        a = int(phenotype[line[0]][horse[header[i2]]]["het"]) + 1
                        phenotype[line[0]][horse[header[i2]]]["het"] = a
                    elif v2 == "1/1" or v2 == "1|1":
                        a = int(phenotype[line[0]][horse[header[i2]]]["hom"]) + 1
                        phenotype[line[0]][horse[header[i2]]]["hom"] = a
                    elif v2 == "./.":
                        a = int(phenotype[line[0]][horse[header[i2]]]["miss"]) + 1
                        phenotype[line[0]][horse[header[i2]]]["miss"] = a

for key in phenotype.keys():
    if phenotype[key] == "NA":
        next
    else:
        for key1 in phenotype[key].keys():
            count_A = 0
            count_R = 0
            for key2 in phenotype[key][key1].keys():
                if key2 == "hom_WT":
                    count_R = count_R + (2*(int(phenotype[key][key1][key2])))
                elif key2 == "het":
                    count_R = count_R + int(phenotype[key][key1][key2])
                    count_A = count_A + int(phenotype[key][key1][key2])
                elif key2 == "hom":
                    count_A = count_A + (2*(int(phenotype[key][key1][key2])))
            phenotype[key][key1]["AC"] = count_A
            phenotype[key][key1]["AF"] = (count_A/(count_A + count_R))

with open("known_variants_summary_by_breed.txt", "w") as outputf:
    print("Phenotype\tBreed\tGenotype\tGenotype_count\tAC\tAF", file = outputf)
    for key in phenotype.keys():
        if phenotype[key] == "NA":
            print(key, "NA\tNA\tNA\tNA\tNA", sep = "\t", file = outputf)
        else:
            breed_AF = 0
            a_count = 0
            for key2 in phenotype[key].keys():
                for key3 in phenotype[key][key2].keys():
                    if key3 == "het" or key3 == "hom":
                        print(key, key2, key3, phenotype[key][key2][key3], phenotype[key][key2]["AC"], phenotype[key][key2]["AF"], sep = "\t", file = outputf)
