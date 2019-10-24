#Common breed rare population
common_breed = {}
with open("common_breed_rare_pop_shared_variants.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[0] + ":" + line[1]
        common_breed[ab] = line[2] + ":" + line[3] + ":" + line[4] + ":" + line[5]
 
with open("common_breed_rare_pop_snpeff_variants.txt","r") as input_file, open("common_breed_rare_pop_snpeff_info.txt", "w") as output_file:
    print("chromosome\tposition\tpop_AC\tpop_AF\tconsequence\timpact\tgene\tcoding\tprotein\tbreed_common\tbreed_AC\tMultiallelic", file = output_file)
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
        if "splice" in consequence:
            consequence = "splice_region_variant"
        elif "5_prime" in consequence:
            consequence = "5_prime_UTR_variant"
        elif "frameshift" in consequence:
            consequence = "frameshift_variant"
        elif "start_lost" in consequence:
            consequence = "start_lost"
        elif "stop_gained" in consequence:
            consequence = "stop_gained"
        elif "stop_lost" in consequence:
            consequence = "stop_lost"
        elif "gene_fusion" in consequence:
            consequence = "gene_fusion_variant"
        elif "synonymous" in consequence:
            consequence = "synonymous_variant"
        elif "missense" in consequence:
             consequence = "missense_variant"
        elif "non_coding" in consequence:
            consequence = "non_coding"
        impact = cd[2]
        gene = cd[3]
        gene = gene.split("-")
        if "CHR_START" in gene[0]:
            gene = gene[-1]
        else:
            if "exon" not in gene[0]:
                gene = gene[0]
            else:
                if "id" in gene[1]:
                    gene = gene[2]
                else:
                    gene = gene[1]
        if cd[9] == "":
            coding = "NA"
        else:
            coding = cd[9]
        if cd[10] == "":
            protein = "NA"
        else:
            protein = cd[10]
        if chrom_pos in common_breed.keys():
            breed_info = common_breed[chrom_pos]
            breed_info = breed_info.split(":")
            print(line[0], line[1], AC, breed_info[2], consequence, impact, gene, coding, protein,breed_info[0],breed_info[1],breed_info[3], file = output_file, sep = "\t")

#Rare breed common population
rare_breed = {}
a = 0
with open("../rare_breed_common_pop/rare_breed_common_pop_shared_variants.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if len(line) > a:
            a = len(line)
            b = line    
        ab = line[0] + ":" + line[1]
        rare_breed[ab] = ":".join(line[2:])

#Figure out what the max length of 
c = list("N"*19)
with open("../rare_breed_common_pop/rare_breed_common_pop_snpeff_variants.txt","r") as input_file, open("../rare_breed_common_pop/rare_breed_common_pop_snpeff_info.txt", "w") as output_file:
    print("chromosome\tposition\tpop_AC\tpop_AF\tconsequence\timpact\tgene\tcoding\tprotein\tbreed_rare\tbreed_AC\tMultiallelic\tother_breed_details\t","\t".join(c),sep="", file = output_file)
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
        if "splice" in consequence:
            consequence = "splice_region_variant"
        elif "5_prime" in consequence:
            consequence = "5_prime_UTR_variant"
        elif "frameshift" in consequence:
            consequence = "frameshift_variant"
        elif "start_lost" in consequence:
            consequence = "start_lost"
        elif "stop_gained" in consequence:
            consequence = "stop_gained"
        elif "stop_lost" in consequence:
            consequence = "stop_lost"
        elif "gene_fusion" in consequence:
            consequence = "gene_fusion_variant"
        elif "synonymous" in consequence:
            consequence = "synonymous_variant"
        elif "missense" in consequence:
             consequence = "missense_variant" 
        elif "non_coding" in consequence:
            consequence = "non_coding"
        impact = cd[2]
        gene = cd[3]
        gene = gene.split("-")
        if "CHR_START" in gene[0]:
            gene = gene[-1]
        else:
            if "exon" not in gene[0]:
                gene = gene[0]
            else:
                if "id" in gene[1]:
                    gene = gene[2]
                else:
                    gene = gene[1]
        if cd[9] == "":
            coding = "NA"
        else:
            coding = cd[9]
        if cd[10] == "":
            protein = "NA"
        else:
            protein = cd[10]
        if chrom_pos in rare_breed.keys():
            breed_info = rare_breed[chrom_pos]
            breed_info = breed_info.split(":")
            a = 24 - len(breed_info)
            b = list("N")
            c = b * a 
            print(line[0], line[1], AC, breed_info[2], consequence, impact, gene, coding, protein,breed_info[0],breed_info[1],"\t".join(breed_info[3:]),"\t".join(c), file = output_file, sep = "\t")
