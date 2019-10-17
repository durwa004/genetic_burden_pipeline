with open("thesis_intersect_annovar.exonic_variant_function", "r") as input_file, open("annovar_variant_type_coding.txt", "w") as output_file:
    print("#CHROM\tPOS\tAC\tAF\tConsequence\tImpact\tGene", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NC_" in line[2]:
            ab = line[17].split("AC=")
            bc = ab[1].split(";")
            AC = bc[0]
            cd = line[1].split(",")
            gene = cd[0]
            if "frameshift" in line[1] or "stopgain" in line[1] or "stoplost" in line[1]:
                impact = "HIGH"
            elif "nonframeshift" in line[1] or "nonsynonymous" in line[1]:
                impact = "MODERATE"
            elif "synonymous SNV" in line[1]:
                impact = "LOW"
            else:
                impact = "UNKNOWN"
            print(line[2], line[3],AC,line[7],line[0],"NA", gene, file = output_file, sep  = "\t")
        else:
            next
