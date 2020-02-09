import os
import gzip

#Get details of frequency of known variants
#Pull out the remapped SNP locations
#Need to get disease for each QTL as well as the tabix code
import numpy as np
directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/known_qtls/split_files/"
chrom_pos = {}
chrom_pos_2 = {}
unmappable = {}
count_m = 0
count_u = 0
count_2 = 0
count = 0
for filename in os.listdir(directory):
    if filename.endswith("_report"):
        with open(directory + "/" + filename) as input_file:
            input_file.readline()
            for line in input_file:
                line = line.rstrip("\n").split("\t")
                count +=1
                if len(line) > 13:
                    if "First Pass" in line[16]:
                        count_m +=1
                        a = line[3] + ":" + line[8]
                        b = line[4] + ":" + line[13]
                        if a in chrom_pos.keys():
                            c = chrom_pos[a] + "," + b
                            chrom_pos[a] = c
                        else:
                            chrom_pos[a] = b
                    elif "Second Pass" in line[16]:
                        count_2 +=1
                        a = line[3] + ":" + line[8]
                        b = line[4] + ":" + line[13]
                        if a in chrom_pos_2.keys():
                            c = chrom_pos_2[a] + "," + b
                            chrom_pos_2[a] = c
                        else:
                            chrom_pos_2[a] = b
                else:
                    unmappable[line[3]] = line[5]
                    count_u +=1

np.setdiff1d(list(chrom_pos_2.keys()), list(chrom_pos.keys()))
chrom_pos["NC_009157:1368081"] = "NC_009172.3:1876173"
chrom_pos["NC_009157:1388861"] = "NC_009172.3:1876173"
chrom_pos["NC_009157:1427118"] = 'NC_009172.3:1876173'


#Need to get disease variants are associated with

#Print out list to extract variants from snpeff file
count = 0
count_2 = 0
with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/QTLs_remapped.sh", "w") as f:
    for filename in os.listdir(directory + "/../"):
        if filename.endswith(".txt"):
            fn = filename.split(".txt")
            with open(directory + "/../" + filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    chrom = line[0].split(".")
                    a = chrom[0] + ":" + line[1]
                    if a in chrom_pos.keys():
                        b = chrom_pos[a].split(",")
                        if len(b) >1:
                            count +=1
                            del chrom_pos[a]
                        else:
                            count_2 +=1
                            ab = chrom_pos[a].split(":")
                            print("/home/mccuem/shared/.local/bin/tabix /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz ", chrom_pos[a], "-", ab[1], " > ", fn[0], ".txt", sep = "", file =f)


#Get dictionary of phenotype and original chrom/pos
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/known_variants/known_QTLs_from_tabix/"

#Get breed info
horse_breed = {}
with open("../../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"
                    
#Get list of horse ids in order of vcf.
header = []
with gzip.open("../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
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

#Get detailed names and symbol
detailed = {} 
with open(path +"/../animalgenomeQTL_for_extraction.txt", encoding = "ISO-8859-1") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1] + ":" + line[2] + ":" + line[5]
        detailed[line[8]] = a

#Get rsid from grep file
rs = {}
with open(path + "/../../dbsnp/EVD_dbsnp/known_qtls/grep_get_rs_chrom_pos.sh", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split()
        rs[line[-1]] = line[1]
#Get variant details for each horse
#with open("known_QTL_locations/No_QTLs_present.txt", "w") as output_file, open("known_QTL_locations/known_QTLs_present.txt", "w") as output2:
with open(path + "/../QTLs_table_tidy.txt", "w") as output2:
#    print("Phenotype", "chrom", "pos", "ref", "alt", "AC", "AF", "\t".join(header[9:]), sep = "\t", file = output2)
#    print("NA", "NA", "NA", "NA", "NA", "NA", "NA", "\t".join(breed), sep = "\t", file = output2)
    print("Phenotype", "Phenotype_abb", "QTL_id", "chrom", "pos", "SNP_EC2", "SNP_EC3","ref", "alt", "AC", "AF", sep = "\t", file = output2)
    for filename in os.listdir(path):
        if filename.endswith(".txt"):
            with open(path + "/" + filename, "r") as input_file:
                if os.stat(path + "/" + filename).st_size !=0:
                    for line in input_file:
                        genotype = []
                        count = 0
                        line = line.rstrip("\n").split("\t")
                        for i in range(len(line)):
                            if "0/1" in line[i] or "0/2" in line[i] or "0/3" in line[i]:
                                genotype.append("1")
                                count +=1
                            elif "1/1" in line[i] or "2/2" in line[i] or "2/1" in line[i] or "1/2" in line[i] or "1/3" in line[i] or "2/3" in line[i] or "3/1" in line[i] or "3/2" in line[i] or "3/3" in line[i]:
                                genotype.append("2")
                                count +=1
                            elif "./." in line[i]:
                                genotype.append("Missing")
                            elif "0/0" in line[i]:
                                genotype.append("0")
                                count +=1
                            elif "/" in line[i]:
                                print(line[i])
                        AC = 0
                        for i in range(len(genotype)):
                            if genotype[i] == "1" or genotype[i] == "2":
                                AC += int(genotype[i])
                        AF = AC/(count*2)
                        rsid = rs[filename]
                        ab = detailed[rsid].split(":")
                        print(ab[2], ab[1], ab[0], line[0], line[1],line[2], rsid, line[3], line[4],AC,AF,sep = "\t", file = output2)
                elif os.stat(path + "/" + filename).st_size == 0:
                    AF = AC/(count*2)
                    rsid = rs[filename]
                    ab = detailed[rsid].split(":")
                    print(ab[2], ab[1], ab[0], ab[3], "NA", "NA", rsid_id, "NA", "NA", "NA", "NA", file = output2, sep = "\t") 

header = []
breed = []
phenotype = {}
count = 0
AF = 0
max_AF = 0
min_AF = 100
with open("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/known_QTLs_present.txt", "r") as input_file, open("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/QTLs_table.txt", "w") as f:
    print("Phenotype\thorse\tbreed\tgenotype", file = f)
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
            AF += float(line[6])
            if float(line[6]) > float(max_AF):
                max_AF = line[6]
            elif float(line[6]) < float(min_AF):
                min_AF =line[6]
            count +=1
            a = line[0].split("_")
            if a[0] in phenotype.keys():
                b = phenotype[a[0]].split(",")
                c = a[0] + "_" + str(len(b))
                d = phenotype[a[0]] + "," + a[0]
                phenotype[a[0]] = d
            else:
                phenotype[a[0]] = a[0]
                c = a[0] + "_0"
            for i in range(len(line)):
                    if line[i] == "1" or line[i] == "2":
                        print(c, header[i], breed[i], line[i], sep = "\t", file = f)

#Mean AF
print(AF/count)
print(max_AF)
print(min_AF)
