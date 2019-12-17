import os
#Get list of crhom:pos to submit to ncbi remapping tool
path = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/known_qtls"

with open("QTL_EC2_chrom_pos.bed", "w") as f:
    for filename in os.listdir(path):
        if "QTL_EC2_chrom_pos" in filename:
            continue
        elif filename.endswith(".txt"):
            with open(filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    print(line[0], line[1], line[1], file = f, sep = "\t")
