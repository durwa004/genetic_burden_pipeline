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
            help="Path to dir with output snpeff files  [required]")
    parser.add_argument(
            "-i", "--ids",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Horse/breed file [required]")
    return parser


if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    horse_ids = args.ids

    horse_breed = {}
    with open(horse_ids, "r") as input_file:
        input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            horse_breed[line[0]] = line[1]
    horse_breed['TWILIGHT'] = "TB"

    breeds = list(set(horse_breed.values()))

    for filename in os.listdir(data):
        if filename.endswith("_common_snpeff.vcf.gz"):
            with gzip.open(data + "/" + filename, "rt") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    header = []
                    if "#CHROM" in line[0]:
                        for i in range(len(line)):
                            header.append(line[i])
                        header1 = header[9:]
                        break

    for item in range(len(breeds)):
        with open(data + "/" + breeds[item] + "_rare_other_breed_common.txt", "w") as output_file:
            print("Common_breed\tCHROM\tPOS\tREF\tALT\tAC\tAF\tconsequence\timpact\tgene\tcoding\tprotein\tlof", "\t".join(header[9:]), sep = "\t", file = output_file)
            for filename in os.listdir(data):
                if filename.endswith("_common_snpeff.vcf.gz"):
                    a = filename.split("_rare_")
                    if a[0] == breeds[item]:
                        with gzip.open(data + "/" + filename, "rt") as input_file:
                            for line in input_file:
                                line = line.rstrip("\n").split("\t")
                                if "#" in line[0]:
                                    next
                                else:
                                    ab = line[7].split(";")
                                    cd = ab[1].split("AF=")
                                    AF = cd[1]
                                    if "," in AF:
                                        ef = AF.split(",")
                                        AF = ef[0]
                                    bc = ab[0].split("AC=")
                                    AC = bc[1]
                                    if "," in AC:
                                        gh = AC.split(",")
                                        AC = gh[0]
                                    de = line[7].split("ANN=")
                                    bc = de[1].split("|")
                                    consequence = bc[1]
                                    coding = bc[9]
                                    protein = bc[10]
                                    impact = bc[2]
                                    gene = bc[3]
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
                                    if "frameshift" in consequence or "start_lost" in consequence or "stop_gained" in consequence or "stop_lost" in consequence:
                                        lof = "y"
                                    else:
                                        lof = "n"
                        print(a[1],line[0], line[1], line[3],line[4], AC, AF, consequence, impact, gene, coding, protein, lof,"\t".join(line[9:]), sep = "\t", file = output_file)
