import gzip
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/"

#Get AF all variants
all_variants = {}
with gzip.open(path + "thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file:
    for line in input_file:
        if "#" in line:
            input_file.readline()
        else:
            line = line.rstrip("\n").split("\t")
            a = line[0] + ":" + line[1]
            b = line[7].split("AF=")
            c = b[1].split(";")
            if "," in c[0]:
                d = c[0].split(",")
                d = d[0]
            else:
                d = c[0]
            if "e-" in d:
                e = d.split("e-")
                f = "0" * int(e[1])
                g = e[0].split(".")
                g = g[0] + g[1]
                d = "0." + f + g
            all_variants[a] = d

#Get AF GB variants
gb_variants = {}
with open(path + "../gb_analysis/high_moderate_variants/genetic_burden_535_horses.txt", "r") as input_file:
    for line in input_file:
        if "#" in line:
            input_file.readline()
        else:
            line = line.rstrip("\n").split("\t")
            a = line[0] + ":" + line[1]
            b = line[7].split("AF=")
            c = b[1].split(";")
            if "," in c[0]:
                d = c[0].split(",")
                d = d[0]
            else:
                d = c[0]
            if "e-" in d:
                e = d.split("e-")
                f = "0" * int(e[1])
                g = e[0].split(".")
                g = g[0] + g[1]
                d = "0." + f + g
            gb_variants[a] = d

#Remove GB variants from all variants
for key,value in gb_variants.items():
    if key in all_variants.keys():
        del all_variants[key]

#Print list of AFs for GB / all variants
with open(path + "../gb_analysis/AF_gb_cf_all_variants/AF_gb_variants.txt", "w") as output_file:
    print("GB_AF", sep = "\t", file = output_file)
    for key,value in gb_variants.items():
        print(value, file = output_file)

with open(path + "../gb_analysis/AF_gb_cf_all_variants/AF_all_variants.txt", "w") as output_file:
    print("all_AF", sep = "\t", file = output_file)
    for key,value in all_variants.items():
        print(value, file = output_file)
