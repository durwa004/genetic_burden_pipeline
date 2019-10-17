import os
import gzip
for filename in os.listdir("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect"):
    if filename.endswith("_AF_less_than_0.005.vcf.gz"):
        with gzip.open(filename, "rt") as input_file:
            for line in input_file:
                line = line.rstrip("\n").split("\t")

