#Common breed rare population
common_breed = {}
with open("common_breed_rare_pop_shared_variants.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[0] + ":" + line[1]
        common_breed[ab] = line[2] + ":" + line[3] + ":" + line[4] + ":" + line[5]
 
with open("common_breed_rare_pop_snpeff_variants.txt","r") as input_file, open("common_breed_rare_pop_snpeff_info.txt", "w") as output_file:
    print("chromosome\tposition\tpop_AC\tpop_AF\tconsequence\timpact\tgene\tcoding\tprotein\tnuclear\tbreed_common\tbreed_AC\tMultiallelic", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        chrom_pos = line[0] +":" + line[1]
        ab = line[7].split(";")
        bc = ab[0].split("AC=")
        if "," in bc[1]:
            cd = bc[1].split(",")
            AC = cd[0]
        else:
            AC = ab[1]
        cd = line[7].split("|")
        consequence = cd[1]
        impact = cd[2]
        gene = cd[3]
        for i in range(len(cd[0:12])):
            if "c." in cd[i]:
                coding = cd[i]
            elif "p." in cd[i]:
                protein = cd[i]
            elif "n." in cd[i]:
                nuclear = cd[i]
            else:
                coding = "NA"
                protein = "NA"
        if chrom_pos in common_breed.keys():
            breed_info = common_breed[chrom_pos]
            breed_info = breed_info.split(":")
            print(line[0], line[1], AC, breed_info[2], consequence, impact, gene, coding, protein,nuclear,breed_info[0],breed_info[1],breed_info[3], file = output_file, sep = "\t")
        else:
            print(chrom_pos)

#Rare breed common population
rare_breed = {}
with open("../rare_breed_common_pop/rare_breed_common_pop_shared_variants.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[0] + ":" + line[1]
        rare_breed[ab] = ":".join(line[2:])

with open("../rare_breed_common_pop/rare_breed_common_pop_snpeff_variants.txt","r") as input_file, open("../rare_breed_common_pop/rare_breed_common_pop_snpeff_info.txt", "w") as output_file:
    print("chromosome\tposition\tpop_AC\tpop_AF\tconsequence\timpact\tgene\tcoding\tprotein\tnuclear\tbreed_rare\tbreed_AC\tMultiallelic\tother_breed_details", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        chrom_pos = line[0] +":" + line[1]
        ab = line[7].split(";")
        bc = ab[0].split("AC=")
        if "," in bc[1]:
            cd = bc[1].split(",")
            AC = cd[0]
        else:
            AC = ab[1]
        cd = line[7].split("|")
        consequence = cd[1]
        impact = cd[2]
        gene = cd[3]
        for i in range(len(cd[0:12])):
            if "c." in cd[i]:
                coding = cd[i]
            elif "p." in cd[i]:
                protein = cd[i]
            elif "n." in cd[i]:
                nuclear = cd[i]
            else:
                coding = "NA"
                protein = "NA"
        if chrom_pos in rare_breed.keys():
            breed_info = rare_breed[chrom_pos]
            breed_info = breed_info.split(":")
            print(line[0], line[1], AC, breed_info[2], consequence, impact, gene, coding, protein,nuclear,breed_info[0],breed_info[1],"\t".join(breed_info[3:]), file = output_file, sep = "\t")
        else:
            print(chrom_pos)
