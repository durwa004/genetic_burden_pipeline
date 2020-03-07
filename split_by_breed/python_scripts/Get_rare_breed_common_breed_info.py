import gzip
import os

directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/"

with open("rare_breed_common_breed.txt", "w") as output_file:
    print("#Breed_rare\tBreed_common\tREF\tALT\tAC_rare\tAC_common\tPop_AF\tMulti-Allelic", file = output_file)
    for entry in os.listdir(directory):
        if os.path.isdir(directory + entry) == True:
            for filename in os.listdir(directory + entry):
                if filename == "0002.vcf.gz":
                    breed = entry.split("_")
                    breed_r = breed[0]
                    breed_c = breed[2]
                    chrom_r = []
                    pos_r = []
                    ref = []
                    alt = []
                    AC_r = []
                    ma = []
                    AF = []
                    with gzip.open(directory + entry + "/" + filename, "rt") as input_file:
                        for line in input_file:
                            line = line.rstrip("\n").split("\t")
                            if "#" in line[0]:
                                next
                            else: 
                                chrom_r.append(line[0])
                                pos_r.append( line[1])
                                ref.append(line[3])
                                alt.append(line[4])
                                a = line[7].split(";")
                                ab = a[0].split("AC=")
                                if "," in ab[1]:
                                    ma.append("Yes")
                                    cd =ab[1].split(",")
                                    AC_r.append(cd[0])
                                else:
                                    ma.append("No")
                                bc = a[1].split("AF=")
                                if "," in bc[1]:
                                    bc = bc[1].split(",")
                                    AF.append(bc[0])
                                else:
                                    AF.append(bc[1])
                elif filename == "0003.vcf.gz":
                    chrom_c = []
                    pos_c = []
                    AC_c = []
                    with gzip.open(directory + entry + "/" + filename, "rt") as input_file:
                        for line in input_file:
                            line = line.rstrip("\n").split("\t")
                            if "#" in line[0]:
                                next
                            else:
                                chrom_c.append(line[0])
                                pos_c.append(line[1])
                                a = line[7].split(";")
                                ab = a[0].split("AC=")
                                if "," in ab[1]:
                                    c = ab[1].split(",")
                                    AC_c.append(c[0])
                                else:
                                    AC_c.append(ab[1])
            for i,v in enumerate(chrom_r):
                for item,value in enumerate(chrom_c):
                    if v == value and pos_r[i] == pos_c[item]:
                        print(breed_r,breed_c,v,pos_r[i],ref[i],alt[i],AC_r[i],AC_c[item],AF[i],ma[i], file = output_file)

path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_chrom_pos/all_shared_variants/"
different = []
with open(path + "/rare_common_breed_chrom_pos_with_breed_details.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        a = line[2].split(":")
        if item == "rare_common_breed_chrom_pos":
            next
        else:
            different.append(len(a))
