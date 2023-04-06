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







################This may be redundant########
with open("No_variants_present.txt", "w") as output_file, open("known_variants_present.txt", "w") as output2:
    print(, "n_hom","AC", "AF", "\t".join(header[9:]), sep = "\t", file = output2)
    for filename in os.listdir("known_variants/"):
        if filename.endswith(".txt"):
            with open("known_variants/" + filename, "r") as input_file:
                if os.stat"known_variants/tabix_files/" + filename).st_size !=0:
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
                elif os.stat(directory + "/tabix_files/" + filename).st_size == 0:
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

#with open(directory + "/../known_disease_locations_2020/known_variants_locations.txt", encoding = "ISO-8859-1") as input_file:
with open(directory + "OMIA_SNVs_05_14_20.txt", encoding = "ISO-8859-1") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[15] == "yes" or line[15] == "y":
            disease["y"][line[11]] = "NA"
        elif line[15] == "no" or line[15] == "n":
            disease["n"][line[11]] = "NA"
        else:
            print(line[15])
        if line[16] == "y":
            causative["y"][line[11]] = "NA"
        elif line[16] == "n":
            causative["n"][line[11]] = "NA"
        else:
            print(line[16])

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
with open(directory + "/known_variants_present_exact_locations_May_2020.txt", "r") as input_file, open(directory + "/variants_by_indvidual.txt", "w") as dz_file:
    print("Phenotype\thorse\tbreed\tgenotype\tAF\tdisease\tcausative", file = dz_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "Phenotype" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
        elif "NA" in line[0]:
            for i in range(len(line)):
                breed.append(line[i])
            next
        else:
            phen = line[0].split(".txt")
            phen = phen[0]
            phenotype.append(line[0])
            if phen in disease["y"].keys():
                dz = "y"
                AF_dz += float(line[8])
                count_dz +=1
                if float(line[8]) > float(max_AF_dz):
                    max_AF_dz = line[8]
                elif float(line[8]) < float(min_AF_dz):
                    min_AF_dz = line[8]
            elif phen in disease["n"].keys():
                dz = "n"
                AF_non_dz += float(line[8])
                count_non_dz +=1
                if float(line[8]) > float(max_AF_non_dz):
                    max_AF_non_dz = line[8]
                elif float(line[8]) < float(min_AF_non_dz):
                    min_AF_non_dz = line[8]
            else:
                print(line)
            if phen in causative["y"].keys():
                cau = "y"
                cau_AF += float(line[8])
                count_cau +=1
                if float(line[8]) > float(max_AF_cau):
                    max_AF_cau = line[8]
                elif float(line[8]) < float(min_AF_cau):
                    min_AF_cau = line[8]
            elif phen in causative["n"].keys():
                cau = "n"
                non_cau_AF += float(line[8])
                count_non_cau += 1
                if float(line[8]) > float(max_AF_non_cau):
                    max_AF_non_cau = line[8]
                elif float(line[8]) < float(min_AF_non_cau):
                    min_AF_non_cau = line[8]
            for i in range(len(line)):
                if line[i] == "het" or line[i] == "hom":
                    print(phen, header[i], breed[i], line[i],line[8],dz,cau, sep = "\t", file = dz_file)
            

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
with open(directory + "/variants_by_indvidual.txt", "r") as input_file,open(directory + "/variants_by_individual_R.txt", "w") as f:
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
