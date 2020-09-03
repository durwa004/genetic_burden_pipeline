import gzip
import os

#Will need to figure out the exact paths for this script
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/lof/"
#path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/"

#Get list of horse ids in order of vcf.
het = {}
hom = {}
header = []
with gzip.open("../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                het[line[i]] = 0
                hom[line[i]] = 0
                header.append(line[i])
            break
header = header[9:]

#Get breed info
horse_breed = {}
with open("../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"
           
variant_caller = {}
#with open(path + "/lof_details.txt", "r") as input_file, open(path + "/lof_by_individual.txt", "w") as output_file:
with open(path + "/genetic_burden_details.txt", "r") as input_file, open(path + "/genetic_burden_by_individual.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[1] + ":" + line[2]
        cd = line[9] + ":" + line[15]
        variant_caller[ab] = cd
        line1 = line[19:]
        for i in range(len(line1)):
             if line1[i] == "het":
                 c = int(het[header[i]]) +1
                 het[header[i]] = c
             elif line1[i] == "hom":
                 c = int(hom[header[i]]) +1
                 hom[header[i]] = c
    for key in het.keys():
        if key in horse_breed.keys():
            print(key, horse_breed[key], het[key], hom[key], sep = "\t", file = output_file)
