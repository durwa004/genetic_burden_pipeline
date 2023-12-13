import argparse
import os
import gzip

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="input vcf.gz file [required]")
    parser.add_argument(
            "-p", "--prog",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Variant annotator: SnpEff or VEP [required]")
    return parser

if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    prog = args.prog
    if prog == "SnpEff":
        prog2 = "ANN="
    elif prog == "VEP":
        prog2 = "CSQ="

#Goal is to convert the snpeff output to a useable text file for analysis, so that we can get a union file of snpeff and annovar output.
with gzip.open(data, "rt") as input_file, open(f"{prog}.hml.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            next
        else:
            details = {}
            if "," in line[4]: 
                chr_pos = line[0] + ":" + line[1] + ":" + line[3] #chr:pos:ref
                details[chr_pos] = {}
                al = line[4].split(",")
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                bc = bc[1].split(",")
                cd = ab[1].split("AF=")
                cd = cd[1].split(",")
                for i in range(len(al)):
                    xv = al[i] + ":" + bc[i] + ":" + cd[i]  #alt:AC:AF
                    details[chr_pos][xv] = {}
                    ef = line[7].split(prog2)
                    ef = ef[1].split(",")
                    for i in range(len(ef)):
                        fg = ef[i].split("|")
                        for key in details[chr_pos].keys():
                            plant = key.split(":")
                            if fg[0] == plant[0]: 
                                plant2 = ":".join(fg) #Allele:annotation:impact:gene_name:gene_ID:Feature:Feature_ID:transcript_biotype:rank/total:HGVSc:HGVSp:cDNA_po/cDNA_ln:CDS_po/CDS_ln:distance_to_feature:errors
                                details[chr_pos][xv][plant2] = "NA"
                for key in details.keys(): #chr:pos:ref
                    for item in details[key].keys():  #alt:AC:AF
                        alpha = key.split(":")
                        beta = item.split(":")
                        complete = []
                        for item2 in details[key][item].keys(): #Allele:annotation:impact:gene_name:gene_ID:Feature:Feature_ID:transcript_biotype:rank/total:HGVSc:HGVSp:cDNA_po/cDNA_ln:CDS_po/CDS_ln:distance_to_feature:errors\
                            item3 = item2.split(":")
                            for i,v in enumerate(item3):
                                complete.append(v)
                        print("\t".join(alpha), "\t".join(beta), "\t".join(complete), sep = "\t", file = output_file)
            else:
                al = line[4].split(",")
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                bc = bc[1]
                cd = ab[1].split("AF=")
                cd = cd[1]
                ef = line[7].split(prog2)
                ef = ef[1].split(",")
                for i in range(len(ef)):
                    fg = ef[i].split("|")
                    plant2 = ":".join(fg) #Allele:annotation:impact:gene_name:gene_ID:Feature:Feature_ID:transcript_biotype:rank/total:HGVSc:HGVSp:cDNA_po/cDNA_ln:CDS_po/CDS_ln:distance_to_feature:errors
                    details[plant2] = "NA"
                complete = []
                for key in details.keys(): 
                    alpha = key.split(":")
                    for i,v in enumerate(alpha):
                        complete.append(v)
                    print(line[0], line[1], line[3], line[4], bc, cd, "\t".join(complete), sep = "\t", file = output_file)

