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
            help="Path to dir with breed thesis intersect files  [required]")
    return parser


if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    variant_dict = {}
    for filename in os.listdir(data):
        if filename.endswith(".vcf.gz"):
             with gzip.open(data + "/" + filename, "rt") as input_file:
                 print(f"Processing: {filename}")
                 breed = filename.split(".")
                 breed = breed[0].split("_")
                 breed = breed[2]
                 for line in input_file:
                     line = line.rstrip("\n").split("\t")
                     if "#" in line[0]:
                        continue
                     else:
                         a = line[0] + ":" + line[1]
                         if a in variant_dict.keys():
                             b = variant_dict[a] + ":" + breed
                             variant_dict[a] = b
                         else:
                             variant_dict[a] = breed

    with open(data + "/overlap_of_variants_by_breed.txt", "w") as output_file:
        for key, value in variant_dict.items():
            a = key.split(":")
            b = value.split(":")
            b = list(set(b))
            print("\t".join(a), "\t".join(b), sep = "\t", file = output_file)
