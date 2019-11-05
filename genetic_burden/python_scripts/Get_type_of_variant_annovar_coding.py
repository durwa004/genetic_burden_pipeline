with open("annovar/thesis_intersect.exonic_variant_function", "r") as input_file, open("annovar/annovar_variant_type_coding.txt", "w") as output_file:
    print("#CHROM\tPOS\tAC\tAF\tConsequence\tImpact\tGene", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "NC_" in line[3]:
            ab = line[18].split("AC=")
            bc = ab[1].split(";")
            AC = bc[0]
            cd = line[2].split(",")
            gene = cd[0]
            if "frameshift" in line[1] or "stopgain" in line[1] or "stoploss" in line[1] or "splice" in line[1]:
                impact = "HIGH"
            elif "nonframeshift" in line[1] or "nonsynonymous" in line[1]:
                impact = "MODERATE"
            elif "synonymous SNV" in line[1]:
                impact = "LOW"
            else:
                impact = "UNKNOWN"
            print(line[3], line[4],AC,line[8],line[1],impact, gene, file = output_file, sep  = "\t")
        else:
            continue
