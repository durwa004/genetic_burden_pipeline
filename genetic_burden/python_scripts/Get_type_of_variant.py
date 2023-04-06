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
with gzip.open(data, "rt") as input_file, open(f"{prog}.hml.txt", "w") as output_file, open(f"{prog}.genes.txt", "w") as output2:
    high = 0
    mod = 0
    low = 0
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            if prog == "VEP":
                if "##INFO=<ID=CSQ" in line[0]:
                    line = line[0].split("Allele|")
                    line = line[1].split("|")
                    print("chrom\tpos\tref\talt\tAC\tAF",  "\t".join(line[0:30]), sep = "\t", file = output_file)
                else:
                    next
            elif prog == "SnpEff":
                if "##INFO=<ID=ANN" in line[0]:
                    line = line[0].split("Allele | ")
                    line = line[1].split(" | ")
                    for i,v in enumerate(line):
                        if "/" in v:
                            v = v.split(" / ")
                            line[i] = v[0] +"_" + v[1]
                    print("chrom\tpos\tref\talt\tAC\tAF", "\t".join(line),"LOF\tgene_name\tgene_ID\tno_transcripts\tpercent_transcripts_affected", sep = "\t", file = output_file)
                else:
                    next
        else:
            chr_pos = line[0] + ":" + line[1] + ":" + line[3] #chr:pos:ref 
            if "," in line[4]: 
                al = line[4].split(",")
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                bc = bc[1].split(",")
                cd = ab[1].split("AF=")
                cd = cd[1].split(",")
                for i in range(len(al)):
                    xv = al[i] + ":" + bc[i] + ":" + cd[i]  #alt:AC:AF
                    ef = line[7].split(prog2)
                    ef = ef[1].split(";")
                    ef = ef[0].split(",")
                    for item1, value1 in enumerate(ef):
                        fg = value1.split("|")
                        if al[i] == fg[0]:
                            plant = []
                            for item, value in enumerate(fg[1:]):
                                if item == 30:
                                    break
                                elif item == 3 or item == 2 or item == 5:
                                    value = value.split(".")
                                    plant.append(value[0])
                                    print(value[0], file = output2)
                                elif value == "":
                                    plant.append("N")
                                else:
                                    plant.append(value)
                            if fg[2] == "HIGH":
                                high +=1
                            elif fg[2] == "MODERATE":
                                mod +=1
                            elif fg[2] == "LOW":
                                low +=1
                            chr_pos1 = chr_pos.split(":")
                            xv = xv.split(":")
                            aa = line[7].split("LOF=(")
                            if len(aa) >1:
                                bb = aa[1].split(")")
                                bb = bb[0].split("|")
                                LOF = "yes"
                                print("\t".join(chr_pos1), "\t".join(xv), "\t".join(plant), LOF, "\t".join(bb), sep = "\t", file = output_file)
                            else:
                                if prog == "SnpEff":
                                    NL = list("N" *4)
                                    LOF = "no"
                                    print("\t".join(chr_pos1), "\t".join(xv), "\t".join(plant), LOF, "\t".join(NL), sep = "\t", file = output_file)
                                elif prog == "VEP":
                                    print("\t".join(chr_pos1), "\t".join(xv), "\t".join(plant), sep = "\t", file = output_file)

                            break
                        else:
                            next
            else:
                al = line[4].split(",")
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                bc = bc[1]
                cd = ab[1].split("AF=")
                cd = cd[1]
                ef = line[7].split(prog2)
                ef = ef[1].split(";")
                ef = ef[0].split(",")
                fg = ef[0].split("|")
                plant = []
                for item, value in enumerate(fg[1:]):
                    if item == 30:
                        break
                    elif item == 3 or item == 2 or item == 5:
                        value = value.split(".")
                        plant.append(value[0])
                        print(value[0], file = output2)
                    elif value == "":
                        plant.append("N")
                    else:
                        plant.append(value)
                if fg[2] == "HIGH":
                    high +=1
                elif fg[2] == "MODERATE":
                    mod +=1
                elif fg[2] == "LOW":
                    low +=1
                chr_pos1 = chr_pos.split(":")
                aa = line[7].split("LOF=(")
                if len(aa) >1:
                    bb = aa[1].split(")")
                    bb = bb[0].split("|")
                    LOF = "yes"
                    print("\t".join(chr_pos1), line[4], bc, cd, "\t".join(plant), LOF, "\t".join(bb), sep = "\t", file = output_file)
                else:
                    if prog == "SnpEff":
                        NL = list("N" *4)
                        LOF = "no"
                        print("\t".join(chr_pos1), line[4], bc, cd, "\t".join(plant), LOF, "\t".join(NL), sep = "\t", file = output_file)
                    elif prog == "VEP":
                        print("\t".join(chr_pos1), line[4], bc, cd, "\t".join(plant), sep = "\t", file = output_file)
    print(f"High impact variants: {high}\nModerate impact variants: {mod}\nLow impact variants: {low}")
