import gzip
#Will need to figure out the exact paths for this script
path = ""
#Get list of horse ids in order of vcf.
with gzip.open("../thesis_intersect.vcf.gz", "rt") as input_file, open("horse_ids_vcf.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                print(line[i], file = output_file)
            break

#Get breed info
horse = {}
with open("../../../../horse_genomes_breeds_tidy.txt", "r") as input_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse[line[0]] = line[1]
horse['TWILIGHT'] = "TB"

#Get  chrom:pos (key) and impact:impact (value)
variant = {}
with open(path + "snpeff_annovar_combined_intersect_high_mod_chrom_pos.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[0] + ":" + line[1]
        b = line[2] + ":" + line[3]
        phenotype[a] = b


#Keep working on this       
with open(path + "known_variants_AFs.txt", "r"
          ) as input_file, open(path + "known_CC_with_breeds.txt", "w"
          ) as output_file, open(path + "known_dz_with_breeds.txt", "w") as output2:
    print("Hom_het\tdisease\tCHROM\tPOS\tID\tREF\tALT\tAC\tAF\tbreed", file = output_file)
    print("Hom_het\tdisease\tCHROM\tPOS\tID\tREF\tALT\tAC\tAF\tbreed", file = output2)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = phenotype[line[1]]
            a = a.split(":")
            if "CC" in a[1] or "Gait" in a[1]:
                if "0/1" in line[i]:
                    print(1,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output_file)
                elif "1/1" in line[i]:
                    print(2,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output_file)
            else:
                if "0/1" in line[i]:
                    print(1,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t",file = output2)
                elif "1/1" in line[i]:
                    print(2,a[1],"\t".join(line[:7]),horse[vcf[i]], sep = "\t", file = output2)

het_C = 0
hom_C = 0  
het_D = 0
hom_D = 0
with open(path + "known_variants_AFs.txt", "r"
          ) as input_file, open(path + "known_CC_with_breeds_hom_het.txt", "w"
          ) as output_file, open(path + "known_dz_with_breeds_hom_het.txt", "w") as output2:
    print("disease\tcount_het\tcount_hom", file = output_file)
    print("disease\tcount_het\tcount_hom", file = output2)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        for i in range(len(line)):
            a = phenotype[line[1]]
            a = a.split(":")
            if "CC" in a[1] or "Gait" in a[1]:
                if "0/1" in line[i]:
                    het_C +=1
                elif "1/1" in line[i]:
                    hom_C +=1
            else:
                if "0/1" in line[i]:
                    het_D+=1
                elif "1/1" in line[i]:
                    hom_D+=1
        print(a[1],het_D,hom_D, sep = "\t", file = output2)
        print(a[1],het_C,hom_C, sep = "\t", file = output_file)

