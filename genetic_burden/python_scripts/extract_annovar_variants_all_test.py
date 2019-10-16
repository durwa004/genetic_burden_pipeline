import os
import argparse


def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-c", "--chrom",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Chromosome to run the analysis on [required]")
    return parser


if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    chromosome = args.chrom

with open("annovar_variant_type_all_NC_009144_3.txt", "w") as output_file, open("../../SnpEff/SnpEff_coding_tidy.txt", "r") as header, open("NC_009144_3_annovar_intersect.variant_function", "r") as input_file:
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

#Simply get the variant consequence
with open("thesis_intersect_annovar.variant_function", "r") as input_file, open("annovar_variant_type_all.txt", "w") as output_file:
    print("#CHROM\tPOS\tAC\tAF\tConsequence\tImpact\tGene", file = output_file)
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        ab = line[17].split("AC=")
        AC = ab[1]
        bc = line[1].split(",")
        gene = bc[0]
        print(line[2], line[3],AC,line[7],line[8],"NA", gene, file = output_file, sep = "\t")
