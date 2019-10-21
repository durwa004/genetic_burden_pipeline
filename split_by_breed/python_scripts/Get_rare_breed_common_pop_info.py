import gzip
import os


#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/breed_pop_differences/"
with open("rare_breed_common_pop.txt", "w") as output_file:
    print("Breed\tchrom\tpos\tAC\tAF\tmultiallelic", file = output_file)
    for filename in os.listdir(directory):
        if filename.endswith("_common_pop.vcf.gz"):
            with gzip.open(directory + "/" +filename, "rt") as input_file:
                a = filename.split("_")
                breed = a[1]
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if "#" in line[0]:
                        next
                    else:
                        ab = line[7].split(";")
                        bc = ab[0].split("AC=")
                        if "," in bc[1]:
                            cd = bc[1].split(",")
                            AC = cd[0]
                            ma = "yes"
                        else:
                            AC = ab[1]
                            ma = "no"
                        de = ab[1].split("AF=")
                        if "," in de[1]:
                            ef = de[1].split(",")
                            AF = ef[0]
                        else:
                            AF = de[1]
                        chrom = line[0]
                        pos = line[1]
                        print(breed,chrom,pos,AC,AF,ma, file = output_file, sep = "\t")

with open("common_breed_rare_pop.txt", "w") as output_file:
    print("Breed\tchrom\tpos\tAC\tAF\tmultiallelic", file = output_file)
    for filename in os.listdir(directory):
        if filename.endswith("_rare_pop.vcf.gz"):
            with gzip.open(directory + "/" +filename, "rt") as input_file:
                a = filename.split("_")
                breed = a[1]
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if "#" in line[0]:
                        next
                    else:
                        ab = line[7].split(";")
                        bc = ab[0].split("AC=")
                        if "," in bc[1]:
                            cd = bc[1].split(",")
                            AC = cd[0]
                            ma = "yes"
                        else:
                            AC = ab[1]
                            ma = "no"
                        de = ab[1].split("AF=")
                        if "," in de[1]:
                            ef = de[1].split(",")
                            AF = ef[0]
                        else:
                            AF = de[1]
                        chrom = line[0]
                        pos = line[1]
                        print(breed,chrom,pos,AC,AF,ma, file = output_file, sep = "\t")
