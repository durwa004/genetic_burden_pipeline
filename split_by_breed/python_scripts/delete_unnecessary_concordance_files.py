import os

directory = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/concordance_gatk_bcftools/"

with open("concordance_gatk_bctools_tidy.sh", "w") as output_file:
    for entry in os.listdir(directory):
        if os.path.isdir(directory + entry) == True:
            for filename in os.listdir(directory + entry):
                if filename.endswith(".vcf"):
                    print("cd ", directory, entry, file = output_file, sep = "")
                    print("rm ", filename, sep= "", file = output_file)
