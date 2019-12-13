import os
import gzip
import argparse

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to dir containing the ibio output files [required]")
    parser.add_argument(
            "-r", "--rare",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Rare allele frequency e.g. 0.05 (5%) [required]")
    parser.add_argument(
            "-c", "--common",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Common allele frequency e.g. 0.10 (10%) [required]")
    return parser


if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    maxAF = args.rare
    minAF = args.common

    for filename in os.listdir(data):
        if "thesis_intersect_without_" in filename and filename.endswith(".vcf.gz"):
            breed = filename.split("thesis_intersect_without_")
            breed = breed[1].split(".vcf.gz")
            breed = breed[0]
            count_rare = 0
            count_common = 0
            count_other = 0
            with gzip.open(data + "/" + filename, "rt") as input_file, open(data + f"/without_{breed}_rare.vcf", "w") as rare, open(data + f"/without_{breed}_common.vcf", "w") as common:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if "#" in line[0]:
                        print("\t".join(line), file = rare)
                        print("\t".join(line), file = common)
                    else:
                        a = line[7].split(";")
                        b = a[0].split("AC=")
                        AC = b[1]
                        if "," in AC:
                            c = AC.split(",")
                            AC = c[0]
                        d = a[2].split("AN=")
                        AN = d[1]
                        if "," in AN:
                            e = AN.split(",")
                            AN = e[1]
                        AF = int(AC)/int(AN)
                        if float(AF) < float(maxAF):
                            print("\t".join(line), file = rare)
                            count_rare +=1
                        elif float(AF) > float(minAF):
                            print("\t".join(line), file = common)
                            count_common +=1
                        else:
                            count_other +=1
            print(f"Without {breed}:Rare-{count_rare}, Common- {count_common}, Other = {count_other}")       
