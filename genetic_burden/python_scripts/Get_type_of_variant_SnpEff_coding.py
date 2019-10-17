import gzip

#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
with gzip.open("thesis_intersect_snpeff.coding.ann.vcf.gz", "rt") as input_file, open("SnpEff_variant_type_coding.txt", "w") as output_file:
    print("#CHROM\tPOS\tAC\tAF\tConsequence\tImpact\tGene", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            next
        else:
            if "NC_" in line[0]:
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                AC = bc[1]
                cd = ab[1].split("AF=")
                AF = cd[1]
                de = line[7].split("ANN=")
                bc = de[1].split("|")
                consequence = bc[1]
                impact = bc[2]
                gene = bc[4]
                print(line[0], line[1], AC, AF, consequence, impact, gene, file = output_file, sep = "\t")
            else:
                next
