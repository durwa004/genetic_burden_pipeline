import os
import gzip

#Get details of frequency of known variants

directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/OMIA_variants"
#Get breed info
horse_breed = {}
with open(directory + "/../../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"
                    
#Get list of horse ids in order of vcf.
header = []
with gzip.open(directory + "/../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
            break

#Get breed in same order as header
breed = []
for i in range(len(header)):
    if header[i] in horse_breed.keys():
        breed.append(horse_breed[header[i]])

#Get variant details for each horse
with open(directory + "/../known_disease_locations_2020/No_variants_present.txt", "w") as output_file, open(directory + "/../known_disease_locations_2020/known_variants_present.txt", "w") as output2:
    print("Phenotype", "chrom", "pos", "ref", "alt", "n_het", "n_hom","AC", "AF", "\t".join(header[9:]), sep = "\t", file = output2)
    print("NA", "NA", "NA", "NA", "NA", "NA", "NA","NA", "NA", "\t".join(breed), sep = "\t", file = output2)
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            with open(directory + "/" + filename, "r") as input_file:
                if os.stat(directory + "/" + filename).st_size !=0:
                    for line in input_file:
                        genotype = []
                        count = 0
                        line = line.rstrip("\n").split("\t")
                        line1 = line[9:]
                        for i in range(len(line1)):
                            a = line1[i].split(":")
                            if "/" in a[0]:
                                b = a[0].split("/")
                            elif "|" in a[0]:
                                b = a[0].split("|")
                            if b[0] == "0" and b[1] == "0":
                                genotype.append("hom_WT")
                            elif b[0] == "." and b[1] == ".":
                                genotype.append("missing")
                            elif b[0] == "0" and int(b[1]) >0:
                                genotype.append("het")
                            elif int(b[0]) >0 and int(b[0]) >0:
                                genotype.append("hom")
                            else:
                                print(line1[i])
                        het = genotype.count("het")
                        hom = genotype.count("hom")
                        WT = genotype.count("hom_WT")
                        AF = (het + (hom*2))/((het + hom + WT)*2)
                        AC = het + (2*hom)
                        print(filename, line[0], line[1],line[3], line[4],het, hom, AC,AF, "\t".join(genotype),sep = "\t", file = output2)
                elif os.stat(directory + "/" + filename).st_size == 0:
                    print(filename, file = output_file)

#Then need to double check that the variants are the exact position that they are supposed to be
#Did this by hand - rename to known_variants_exact_locations.txt

#Get dictionary of causative vs associated and disease or not.
causative = {}
causative["y"] = {}
causative["n"] = {}
disease = {}
disease["y"] = {}
disease["n"] = {}

with open(directory + "/../known_disease_locations_2020/known_variants_locations.txt", encoding = "ISO-8859-1") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[2] == "y":
            disease["y"][line[1]] = "NA"
        elif line[2] == "n":
            disease["n"][line[1]] = "NA"
        if "ssociated" in line[3]:
            causative["n"][line[1]] = "NA"
        elif "ausative" in line[3]:
            causative["y"][line[1]] = "NA"

print(len(disease["y"].keys()))
print(len(disease["n"].keys()))
print(len(causative["y"].keys()))
print(len(causative["n"].keys()))
print(len(set(disease["y"].keys() & causative["y"].keys())))
print(len(set(disease["y"].keys() & causative["n"].keys()))) 
print(len(set(disease["n"].keys() & causative["y"].keys()))) 
print(len(set(disease["n"].keys() & causative["n"].keys()))) 

header = []
breed = []
phenotype = []
AF_dz = 0
max_AF_dz = 0
min_AF_dz = 100
count_dz = 0
count_non_dz = 0
AF_non_dz = 0
max_AF_non_dz = 0
min_AF_non_dz = 100
count_cau = 0
count_non_cau = 0
min_AF_cau = 0
max_AF_cau = 0
min_AF_non_cau = 0
max_AF_non_cau = 0
cau_AF = 0
non_cau_AF = 0
with open(directory + "/../known_disease_locations_2020/known_variants_present_exact_locations.txt", "r") as input_file, open(directory + "/../known_disease_locations_2020/variants_by_indvidual.txt", "w") as dz_file:
    print("Phenotype\thorse\tbreed\tgenotype\tAF\tdisease\tcausative", file = dz_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "Phenotype" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
            next
        elif "NA" in line[0]:
            for i in range(len(line)):
                breed.append(line[i])
            next
        else:
            phenotype.append(line[0])
            if line[0] in disease["y"].keys():
                dz = "y"
                AF_dz += float(line[8])
                count_dz +=1
                if float(line[8]) > float(max_AF_dz):
                    max_AF_dz = line[8]
                elif float(line[8]) < float(min_AF_dz):
                    min_AF_dz = line[8]
            elif line[0] in disease["n"].keys():
                dz = "n"
                AF_non_dz += float(line[8])
                count_non_dz +=1
                if float(line[8]) > float(max_AF_non_dz):
                    max_AF_non_dz = line[8]
                elif float(line[8]) < float(min_AF_non_dz):
                    min_AF_non_dz = line[8]
            else:
                print(line[0])
            if line[0] in causative["y"].keys():
                cau = "y"
                cau_AF += float(line[8])
                count_cau +=1
                if float(line[8]) > float(max_AF_cau):
                    max_AF_cau = line[8]
                elif float(line[8]) < float(min_AF_cau):
                    min_AF_cau = line[8]
            elif line[0] in causative["n"].keys():
                cau = "n"
                non_cau_AF += float(line[8])
                count_non_cau += 1
                if float(line[8]) > float(max_AF_non_cau):
                    max_AF_non_cau = line[8]
                elif float(line[8]) < float(min_AF_non_cau):
                    min_AF_non_cau = line[8]
            for i in range(len(line)):
                if line[i] == "het" or line[i] == "hom":
                    print(line[0], header[i], breed[i], line[i],line[8],dz,cau, sep = "\t", file = dz_file)
            

print(AF_dz/count_dz)
print(max_AF_dz)
print(min_AF_dz)
print(AF_non_dz/count_non_dz)
print(max_AF_non_dz)
print(min_AF_non_dz)
print(cau_AF/count_cau)
print(max_AF_cau)
print(min_AF_cau)
print(non_cau_AF/count_non_cau)
print(max_AF_non_cau)
print(min_AF_non_cau)

#Convert variants_bt_indvidual.txt to a version for R disease/breed/genotype/count
dz_details = {}
with open(directory + "/../known_disease_locations_2020/variants_bt_indvidual.txt", "r") as input_file,open(directory + "/../known_disease_locations_2020/variants_by_individual_R.txt", "w") as f:
    input_file.readline()
    print("Phenotype\tbreed\tgenotype\tgenotype_count\tAC\tAF\tdisease\tcausative", file = f)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[4] + ":" + line[5] + ":" + line[6]
        if a in dz_details.keys():
            if line[2] in dz_details[a].keys():
                if line[3] == "het":
                    dz_details[a][line[2]]["het"] +=1
                else:
                    dz_details[a][line[2]]["hom"] +=1
            else:
                dz_details[a][line[2]] = {}
                if line[3] == "het":
                    dz_details[a][line[2]]["het"] = 1
                    dz_details[a][line[2]]["hom"] = 0  
                else:
                    dz_details[a][line[2]]["het"] = 0
                    dz_details[a][line[2]]["hom"] = 1
        else:
            dz_details[a] = {}
            dz_details[a][line[2]] = {}
            if line[3] == "het":
                dz_details[a][line[2]]["het"] = 1
                dz_details[a][line[2]]["hom"] = 0
            else:
                dz_details[a][line[2]]["het"] = 0
                dz_details[a][line[2]]["hom"] = 1
    for phenotype in dz_details.keys():
        phenotype1 = phenotype.split(":")
        for breed in dz_details[phenotype].keys():
            total = (dz_details[phenotype][breed]["het"] + 2*(dz_details[phenotype][breed]["hom"]))
            print(phenotype1[0], breed, "het", dz_details[phenotype][breed]["het"],total, phenotype1[1], phenotype1[2], phenotype1[3], sep = "\t", file = f)
            print(phenotype1[0], breed, "hom", dz_details[phenotype][breed]["hom"],total, phenotype1[1], phenotype1[2], phenotype1[3], sep = "\t", file = f)
