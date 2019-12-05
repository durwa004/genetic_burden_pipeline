import gzip
import os

#Will need to figure out the exact paths for this script
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/gb_analysis/high_moderate_variants/"
#Get list of horse ids in order of vcf.
het = {}
hom = {}
missing = {}
header = []
with gzip.open("../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
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
with open("../../../horse_genomes_breeds_all.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"
            
#Get genetic burden per individual (overall)
with open(path + "/genetic_burden_535_horses.txt", "r") as input_file, open(path + "/ann_se_gb_by_individual_breed.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        line1 = line[9:]
        for i in range(len(line1)):
            if "0/1" in line1[i] or "0/2" in line1[i] or "0/3" in line1[i]:
                a = het[header[i]]
                b = int(a) + 1
                het[header[i]] = b
            elif "1/1" in line1[i] or "2/2" in line1[i] or "1/2" in line1[i] or "1/3" in line1[i] or "2/3" in line1[i] or "2/1" in line1[i]:
                a = hom[header[i]]
                b = int(a) + 1
                hom[header[i]] = b
            elif "./." in line1[i]:
                a = missing[header[i]]
                b = int(a) + 1
                missing[header[i]] = b
    for key in het.keys():
        if key in horse_breed.keys():
            print(key, horse_breed[key], het[key], hom[key], missing[key],sep = "\t", file = output_file)

#Get number of lof variants per individual
header1 = ["NA", "NA"] + header
with open(path + "/lof_variants_breed.txt", "r") as input_file, open(path + "/ann_se_lof_by_individual_breed.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            if "0/1" in line[i] or "0/2" in line[i] or "0/3" in line[i]:
                a = het[header1[i]]
                b = int(a) + 1
                het[header1[i]] = b
            elif "1/1" in line[i] or "2/2" in line[i] or "1/2" in line[i] or "1/3" in line[i] or "2/3" in line[i] or "2/1" in line[i]:
                a = hom[header1[i]]
                b = int(a) + 1
                hom[header1[i]] = b
            elif "./." in line[i]:
                a = missing[header1[i]]
                b = int(a) + 1
                missing[header1[i]] = b
    for key in het.keys():
        if key in horse_breed.keys():
            print(key, horse_breed[key], het[key], hom[key], missing[key],sep = "\t", file = output_file)
