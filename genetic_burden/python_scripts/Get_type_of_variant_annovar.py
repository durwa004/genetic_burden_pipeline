import os

with open("../annovar/annovar_coding.txt", "r") as input_file, open("SnpEff_coding_tidy.txt", "r") as header, open("annovar_variant_type_all.txt", "w") as output_file:
    info = header.readline()
    print(info, file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        consequence = line[1]
        if "frameshift" in line[1] or "stopgain" in line[1] or "stoplost" in line[1]:
            impact = "HIGH"
        elif "nonframeshift" in line[1] or "nonsynonymous" in line[1]:
            impact = "MODERATE"
        elif "synonymous SNV" in line[1]:
            impact = "LOW"
        elif "unknown" in line[1]:
            impact = "UNKNOWN"
        if line[2] == "UNKNOWN":
            b = line[18].split("AC=")
            c = b[1].split(";")
            AC = c[0]
            if "AF_" in line[18]:
                e = line[18].split("AF=")
                f = e[1].split(";")
                AF = f[0]
            else:
                AF = "NA"
            print("\t".join(line[11:18]), AC, AF, consequence, impact, "NA", "NA", "NA", "NA", "\t".join(line[20:]), file = output_file, sep = "\t")
        else:
            a = line[2].split(":")
            if a[2] == "wholegene,":
                gene = a[0]
                coding = "NA"
                protein = "NA"
                b = line[18].split("AC=")
                c = b[1].split(";")
                AC = c[0]
                if "AF_" in line[18]:
                    e = line[18].split("AF=")
                    f = e[1].split(";")
                    AF = f[0]
                else:
                    AF = "NA"
                print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein, "\t".join(line[20:]), file = output_file, sep = "\t")
            else:
                gene = a[0]
                coding = a[3]
                protein = a[4]
                b = line[18].split("AC=")
                c = b[1].split(";")
                AC = c[0]
                if "AF_" in line[18]:
                    e = line[18].split("AF=")
                    f = e[1].split(";")
                    AF = f[0]
                else:
                    AF = "NA"
                print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein, "\t".join(line[20:]), file = output_file, sep = "\t")
