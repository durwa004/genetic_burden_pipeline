import os

directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/"

with open("breed_common_rare_variants_tidy.sh", "w") as output_file:
    for entry in os.listdir(directory):
        if os.path.isdir(directory + entry) == True:
            for filename in os.listdir(directory + entry):
                if filename == "0001.vcf":
                    print("cd ", directory, entry, file = output_file, sep = "")
                    print("rm ", filename, sep= "", file = output_file)
                elif filename == "0000.vcf" or filename == "0003.vcf" or filename == "0002.vcf":
                    print("/home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip ",filename, sep = "", file = output_file)
                    print("/home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix ", filename,".gz", file = output_file, sep = "") 
