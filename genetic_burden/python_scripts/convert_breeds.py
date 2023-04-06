#Convert ibio breeds (from Samples spreadsheet) to breeds for GB analysis

#Get horse IDs
```
module load bcftools
bcftools query
bcftools query -l /panfs/jay/groups/6/durwa004/shared/PopulationVCF/joint_genotype.goldenPath.vep.vcf.gz > horse_IDs.txt
```

#Add in breed info from excel sheet
```
scp durwa004@mesabi.msi.umn.edu://home/durwa004/durwa004/genetic_burden/horse_IDs.txt
#Input IDs from M# database
scp horse_genomes_breeds_all.txt durwa004@mesabi.msi.umn.edu://home/durwa004/durwa004/genetic_burden/
module load python/3.6.3 
ipython
```

with open("horse_genomes_breeds_all.txt", "r") as input_file, open("horse_genomes_breeds_tidy.txt","w") as output_file:
    print("horse_id", "breed",  file = output_file, sep = "\t")
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        horse = line[0]
        breed = line[1]
        if breed == "Arabian" or breed == "Arab":
            print(horse, "Arabian", file = output_file, sep = "\t")
        elif breed == "Belgian":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "Clydesdale":
            print(horse, breed, file = output_file, sep = "\t") 
        elif breed == "Icelandic" or breed == "Icelandic Horse":
            print(horse, "Icelandic", file = output_file, sep = "\t")
        elif breed == "Morgan":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "Quarter Horse" or breed == "QH" or breed == "QHx" or breed == "QH (App)" or breed == "Paint":
            print(horse, "QH", file = output_file, sep = "\t")
        elif breed == "Shetland":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "STB" or breed == "Standardbred" or breed == "StandardBred" or breed == "StandardBred (TROTTER)" or breed == "Standardbred (pacer)" or breed == "Standardbred (trotter)" or breed == "Standardbred (PACER)":
            print(horse, "STB", file = output_file, sep = "\t")
        elif breed == "Thoroughbred" or breed == "TB" or breed == "Thoroughbred x":
            print(horse, "TB", file = output_file, sep = "\t")
        elif breed == "Welsh Pony" or breed == "WP":
            print(horse, "WP", file = output_file, sep = "\t")
        elif breed == "Franchese Montagne":
            print(horse, "FM", file = output_file, sep = "\t")
        elif breed == "Przewalski":
            next
        else:
            print(horse, "Other", file = output_file, sep = "\t")
