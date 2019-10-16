import gzip
import os

#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
with gzip.open("thesis_intersect_snpeff.ann.vcf.gz", "rt") as input_file, open("SnpEff_variant_type_all.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line:
            print("\t".join(line[0:7]), "AC", "AF", "Consequence","Impact", "Gene", "Where", "coding_change", "protein_change", "\t".join(line[9:]), file = output_file, sep = "\t")
        elif "#" in line[0]:
            next
        else:
            ab = line[7].split("ANN=")
            bc = ab[1].split("|")
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
            print("\t".join(line[0:7]), AC, AF, consequence, impact, gene, where, coding, protein, "\t".join(line[9:]), file = output_file, sep = "\t")
