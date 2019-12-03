import gzip

path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/"
#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
with gzip.open(path + "thesis_intersect_af_over_0.5.vcf.gz", "rt") as input_file, open(path + "thesis_intersect_variants_present_in_all.vcf", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            print("\t".join(line), file = output_file)
        else:
            line1 = line[9:]
            count = 0
            for i in range(len(line1)):
                if "0/0" in line1[i]:
                    count +=1
                else:
                    continue
            if count > 0:
                next
            else:
                print("\t".join(line), file = output_file)
