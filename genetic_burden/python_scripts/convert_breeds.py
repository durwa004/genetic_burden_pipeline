#Convert ibio breeds (from Samples spreadsheet) to breeds for GB analysis

with open("horse_genomes_breeds.txt", "r") as input_file, open("horse_genomes_breeds_tidy.txt","w") as output_file:
    print("horse_id", "breed",  file = output_file, sep = "\t")
    input_file.readline()
    for line in input_file:
        horse,breed = line.rstrip("\n").split("\t")
        if breed == "Arabian":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "Belgian":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "Clydesdale":
            print(horse, breed, file = output_file, sep = "\t") 
        elif breed == "Icelandic" or breed == "Icelandic Horse":
            print(horse, "Icelandic", file = output_file, sep = "\t")
        elif breed == "Morgan":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "Quarter Horse" or breed == "QH":
            print(horse, "QH", file = output_file, sep = "\t")
        elif breed == "Shetland":
            print(horse, breed, file = output_file, sep = "\t")
        elif breed == "Standardbred" or breed == "StandardBred" or breed == "StandardBred (TROTTER)" or breed == "Standardbred (pacer)" or breed == "Standardbred (trotter)":
            print(horse, "STB", file = output_file, sep = "\t")
        elif breed == "Thoroughbred":
            print(horse, "TB", file = output_file, sep = "\t")
        elif breed == "Welsh Pony":
            print(horse, "WP", file = output_file, sep = "\t")
        elif breed == "Przewalski":
            next
        else:
            print(horse, "Other", file = output_file, sep = "\t")

#Print out list of ids for bcftools view
with open("horse_genomes_breeds_tidy.txt","r") as input_file, open("joint_gvcf/joint_intersect/horse_ids_all.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        horse,breed = line.rstrip("\n").split("\t")
        print(horse, file = output_file)

