import gzip
import os

#Will need to figure out the exact paths for this script
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/"
#Get list of horse ids in order of vcf.
het = {}
hom = {}
missing = {}
header = []
with gzip.open("../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                het[line[i]] = 0
                hom[line[i]] = 0
                missing[line[i]] = 0
                header.append(line[i])
            break

#Get breed info
horse_breed = {}
with open("../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"
            
#Get genetic burden per individual (overall)
with open(path + "/lof_snpeff.txt", "r") as input_file, open(path + "/lof_by_individual.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = line1[i].split(":")
            if "/" in a[0]:
                b = a[0].split("/")
            elif "|" in a[0]:
                b = a[0].split("|")
            if b[0] == ".":
                miss +=1
            elif b[0] == "0":
                if int(b[1]) > 0:
                c = het[header[i]]
                d = int(c) + 1
                het[header[i]] = d                   
            elif int(b[0]) >0:
                c = hom[header[i]]
                d = int(c) + 1
                hom[header[i]] = d
    for key in het.keys():
        if key in horse_breed.keys():
            print(key, horse_breed[key], het[key], hom[key], missing[key],sep = "\t", file = output_file)
