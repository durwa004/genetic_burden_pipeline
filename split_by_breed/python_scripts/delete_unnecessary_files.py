import os

directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/rare_common_breed_pop/breed_pop_rare_common_vcfs/"

with open("breed_pop_common_rare_tidy.sh", "w") as output_file:
    for entry in os.listdir(directory):
        if os.path.isdir(directory + entry) == True:
            print("cd " + directory + "/" + entry + "/",sep= "", file = output_file)
            for filename in os.listdir(directory + entry):
                if filename.endswith(".vcf"):
                    print("bcftools stats ", filename, " > ", filename, ".stats", file = output_file, sep = "")
                if filename == "0001.vcf": #or filename == "0003.vcf":
                    print("rm ", filename, sep= "", file = output_file)
#Need to include 0003.vcf for the common analysis to get AF
                elif filename == "0000.vcf" or  filename == "0002.vcf" or filename == "0003.vcf":
                    print("/home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip ",filename, sep = "", file = output_file)
                    print("/home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix ", filename,".gz", file = output_file, sep = "") 
