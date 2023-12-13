import gzip
import os

#Get list of horse ids in order of vcf.
vcf = "joint_genotype.goldenPath.vep.subset.snpeff.hml.vcf.gz"
breeds_file = "horse_genomes_breeds_tidy.txt"
header = []
with gzip.open(vcf, "rt") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                header.append(line[i])
            break
header = header[9:]

#Get breed info
horse_breed = {}
with open(breeds_file, "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse_breed[line[0]] = line[1]

#Get variant location           
variants = {}
#with open("SnpEff.VEP.mod.intersect.txt", "r") as input_f:
with open("SnpEff.VEP.intersect.txt", "r") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        a = line[1] + ":" + line[2] + ":" + line[4]
        if a in variants.keys():
            print(line)
        else:
            variants[a] = line[3]

RefSeq = {}
with open("RefSeq_chromosomes.txt", "r") as input_f:
    input_f.readline()
    for line in input_f:
        line = line.rstrip("\n").split("\t")
        RefSeq[line[1]] = line[0]

#Get genotypes
all_gt = []
#with gzip.open(vcf, "rt") as input_file, open("ibio_SnpEff.VEP.mod.intersect.individual.txt", "w") as output_f, open("mod_non_intersect.txt", "w") as outpu2:
with gzip.open(vcf, "rt") as input_file, open("SnpEff.VEP.intersect.individual.txt", "w") as output_f, open("non_intersect.txt", "w") as outpu2:
#with gzip.open(vcf, "rt") as input_file, open("ss.SnpEff.VEP.intersect.individual.txt", "w") as output_f, open("ss_non_intersect.txt", "w") as outpu2:
#with gzip.open(vcf, "rt") as input_file, open("tu.SnpEff.VEP.intersect.individual.txt", "w") as output_f, open("tu_non_intersect.txt", "w") as outpu2:
#with gzip.open(vcf, "rt") as input_file, open("ti.SnpEff.VEP.intersect.individual.txt", "w") as output_f, open("ti_non_intersect.txt", "w") as outpu2:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            header = []
            for i,v in enumerate(line[9:]):
                header.append(v)
            print("CHROM", line[1], line[3], line[4], "\t".join(header), sep = "\t", file = output_f)
#            print("\t".join(header), sep = "\t", file = output_f)
        elif "#" in line[0] or "NW_" in line[0]:
            next
        else:
            if "," in line[4]:
                genotypes = []
                alt = line[4].split(",")
                loc1 = line[0] + ":" + line[1] + ":" + alt[0]
                loc2 = line[0] + ":" + line[1] + ":" + alt[1]
#                loc1 = RefSeq[line[0]] + ":" + line[1] + ":" + alt[0]
#                loc2 = RefSeq[line[0]] + ":" + line[1] + ":" + alt[1]
                if loc1 in variants.keys():
                    for i,v in enumerate(line[9:]):
                        gt = v.split(":")
                        all_gt.append(gt[0])
                        genotypes.append(gt[0])
                    for i in range(len(alt)):
#                        print(RefSeq[line[0]], line[1], line[3], alt[i], "\t".join(genotypes), sep = "\t", file = output_f)
                        print(line[0], line[1], line[3], alt[i], "\t".join(genotypes), sep = "\t", file = output_f)
#                        print("\t".join(genotypes), sep = "\t", file = output_f)
                else:
                    a = line[7].split("AF=")
                    a = a[1].split(";")
                    a = a[0].split(",")
                    for i in range(len(a)):
                        print(a[i], file = outpu2)
            else:
                genotypes = []
#                loc1 = RefSeq[line[0]] + ":" + line[1] + ":" + line[4]
                loc1 = line[0] + ":" + line[1] + ":" + line[4]
                if loc1 in variants.keys():
                    for i,v in enumerate(line[9:]):
                        gt = v.split(":")
                        all_gt.append(gt[0])
                        genotypes.append(gt[0])
#                    print(RefSeq[line[0]], line[1], line[3], line[4], "\t".join(genotypes), sep = "\t", file = output_f)
                    print(line[0], line[1], line[3], line[4], "\t".join(genotypes), sep = "\t", file = output_f)
#                    print("\t".join(genotypes), sep = "\t", file = output_f)
                else:
                    a = line[7].split("AF=")
                    a = a[1].split(";")
                    a = a[0].split(",")
                    print(a[0], file = outpu2)

#Get all possible genotype values for R analysis
print(set(all_gt))






###################OLD########

#find variants unique to breeds  
unique = []
u_b = []
for item in gb.keys():
    a = gb[item].split(":")
    if len(set(a)) <3:
        unique.append(item)
        u_b.append(a[1])

with open(path + "/no_homozygotes_details.txt", "r") as input_file, open(path + "/unique_no_homozygotes.txt", "w") as output_file:
    print("Breed\tCHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = output_file)
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        for i in range(len(unique)):
            if a == unique[i]:
                print(u_b[i], "\t".join(line), file = output_file, sep ="\t")


