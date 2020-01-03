import gzip
import os

#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
with gzip.open("../SnpEff/thesis_intersect_snpeff.coding.ann.vcf.gz", "rt") as input_file, open("../SnpEff/SnpEff_coding_tidy.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line:
            print("\t".join(line[0:7]), "AC", "AF", "Consequence","Impact", "Gene", "Where", "coding_change", "protein_change", "lof", "\t".join(line[9:]), file = output_file, sep = "\t")
        elif "#" in line[0]:
            next
        else:
            lof = "NA"
            ab = line[7].split("ANN=")
            bc = ab[1].split("|")
            ef = line[7].split(";")
            for i in range(len(ef)):
                fe = ef[i].split("LOF=")
                if len(fe) > 1:
                    lof = "y"
                else:
                    lof = "n"
            if bc[0] == line[4]:
                a = line[7].split("AC=")
                b = a[1].split(";")
                AC = b[0]
                if "AF_" in line[7]:
                    c = line[7].split("AF=")
                    d = c[1].split(";")
                    AF = d[0]
                else:
                    AF = "NA"
                consequence = bc[1]
                impact = bc[2]
                gene = bc[3]
                where = bc[7]
                coding = bc[9]
                protein = bc[10]
            else:
                a = line[7].split("AC=")
                b = a[1].split(";")
                cd = b[0].split(",")
                AC =cd[0]
                if "AF_" in line[7]:
                    c = line[7].split("AF=")
                    d = c[1].split(";")
                    de = d[0].split(",")
                    AF = de[0]
                else:
                    AF = "NA"
                consequence = bc[1]
                impact = bc[2]
                gene = bc[3]
                where = bc[7]
                coding = bc[9]
                protein = bc[10]
                ef = line[7].split(";")
            print("\t".join(line[0:7]), AC, AF, consequence, impact, gene, where, coding, protein,lof, "\t".join(line[9:]), file = output_file, sep = "\t")

#Extract annovar variants into same format as snpeff variants
with open("../annovar/annovar_exonic_variant_function/thesis_intersect.exonic_variant_function", "r") as input_file, open("../SnpEff/SnpEff_coding_tidy.txt", "r") as header, open("../annovar/annovar_coding_tidy.txt", "w") as output_file, open("../annovar/annovar_coding_tidy_lof.txt", "w") as output:
    info = header.readline()
    print(info, file = output_file)
    print(info, file = output)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        consequence = line[1]
        de = line[1].split (" ")
        if "frameshift" ==de[0] or "stopgain" in line[1] or "stoploss" in line[1] or "splice" in line[1]:
            impact = "HIGH"
        elif "nonframeshift" in line[1] or "nonsynonymous" in line[1]:
            impact = "MODERATE"
        elif "synonymous SNV" in line[1]:
            impact = "LOW"
        elif "unknown" in line[1]:
            impact = "UNKNOWN"
        else:
            print(line[1])
        if "frameshift" == de[0] or "stopgain" == line[1] or "splice" in line[1]:
            lof = "y"
        else:  
            lof = "n"
        AF = line[8] 
        ab = line[18].split("AC=")
        bc = ab[1].split(";")
        AC = bc[0]
        cd = line[2].split(",")
        gene = cd[0]
        if line[2] == "UNKNOWN":
            print("\t".join(line[11:18]), AC, AF, consequence, impact, "NA", "NA", "NA", "NA", lof, "\t".join(line[20:]), file = output_file, sep = "\t")
            if lof == "y":
                print("\t".join(line[11:18]), AC, AF, consequence, impact, "NA", "NA", "NA", "NA", lof, "\t".join(line[20:]), file = output, sep = "\t")
        else:
            a = line[2].split(":")
            if a[2] == "wholegene,":
                gene = a[0]
                coding = "NA"
                protein = "NA"
                print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein,lof, "\t".join(line[20:]), file = output_file, sep = "\t")
                if lof == "y":
                    print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein,lof, "\t".join(line[20:]), file = output, sep = "\t")
            else:
                gene = a[0]
                coding = a[3]
                protein = a[4]
                print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein, lof, "\t".join(line[20:]), file = output_file, sep = "\t")
                if lof == "y":
                    print("\t".join(line[11:18]), AC, AF, consequence, impact, gene, "NA", coding, protein, lof, "\t".join(line[20:]), file = output, sep = "\t")
