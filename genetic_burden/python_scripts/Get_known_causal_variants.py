import os
import gzip

#Get details of frequency of known variants

#Get breed info
horse_breed = {}
with open("../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]
horse_breed['TWILIGHT'] = "TB"
                    
#Get list of horse ids in order of vcf.
header = []
with gzip.open("../SnpEff/thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
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
with open("known_disease_locations/No_variants_present.txt", "w") as output_file, open("known_disease_locations/known_variants_present.txt", "w") as output2:
    print("Phenotype", "chrom", "pos", "ref", "alt", "AC", "AF", "\t".join(header[9:]), sep = "\t", file = output2)
    print("NA", "NA", "NA", "NA", "NA", "NA", "NA", "\t".join(breed), sep = "\t", file = output2)
    for filename in os.listdir():
        if filename.endswith(".txt"):
            with open(filename, "r") as input_file:
                if os.stat(filename).st_size !=0:
                    for line in input_file:
                        genotype = []
                        count = 0
                        line = line.rstrip("\n").split("\t")
                        for i in range(len(line)):
                            if "0/1" in line[i]:
                                genotype.append("1")
                                count +=1
                            elif "1/1" in line[i]:
                                genotype.append("2")
                                count +=1
                            elif "./." in line[i]:
                                genotype.append("Missing")
                            elif "0/0" in line[i]:
                                genotype.append("0")
                                count +=1
                        AC = 0
                        for i in range(len(genotype)):
                            if genotype[i] == "1" or genotype[i] == "2":
                                AC += int(genotype[i])
                        AF = AC/(count*2)
                        print(filename, line[0], line[1],line[3], line[4],AC,AF, "\t".join(genotype),sep = "\t", file = output2)
                elif os.stat(filename).st_size == 0:
                    print(filename, file = output_file)

#Then need to double check that the variants are the exact position that they are supposed to be
#Did this by hand - rename to known_variants_for_analysis.txt
header = []
breed = []
with open("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/known_variants_for_analysis.txt", 
          "r") as input_file, open("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/dz_variants_table.txt", 
             "w") as dz_file, open("/Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/known_causal_variants/non-dz_variants_table.txt",
                "w") as non_dz_file:
    print("Disease\thorse\tbreed\tgenotype", file = dz_file)
    print("Disease\thorse\tbreed\tgenotype", file = non_dz_file)
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
            if line[2] == "n":
                for i in range(len(line)):
                    if i >8:
                        if line[i] == "1" or line[i] == "2":
                            print(line[1], header[i], breed[i], line[i], sep = "\t", file = non_dz_file)
            else:
                for i in range(len(line)):
                    if i >8:
                        if line[i] == "1" or line[i] == "2":
                            print(line[1], header[i], breed[i], line[i], sep = "\t", file = dz_file)
            
