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
            help="Path to dir with output snpeff file  [required]")
    parser.add_argument(
            "-l", "--lists",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="directory containing list of chrom/pos for each rare/common breed pair [required]")
    return parser


if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    lists = os.path.abspath(args.lists)

    cp_dict = {}
    for filename in os.listdir(lists):
        if os.path.isdir(lists + "/" + filename):
            continue
        else: 
            with open(lists + "/" + filename) as input_file:
                a = filename.split(".txt")
                for line in input_file:
                    chrom, pos = line.rstrip("\n").split("\t")
                    b = chrom + ":" + pos
                    if b in cp_dict.keys():
                        c = cp_dict[b] + "," + a[0]
                        cp_dict[b] = c
                    else:
                        cp_dict[b] = a[0]

#Get chrom/pos for all of the variants that are differentiated - can't get all the details, otherwise there is too much information
    n_shared = []
    genes = []
    AF_list = []
    lof = []
    with gzip.open(data + "/rare_common_breed_coding.vcf.gz", "rt") as input_file, open(data + "/breed_rare_other_breed_common.txt", "w") as f:
        print("Gene\tlof", sep = "\t", file = f)
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            if "#" in line[0]:
                continue 
            else:
                a = line[0] + ":" + line[1]
                if a in cp_dict.keys():
                    b = cp_dict[a].split(",")
                    n_shared.append(len(b))
                ab = line[7].split(";")
                cd = ab[1].split("AF=")
                AF = cd[1]
                if "," in AF:
                    ef = AF.split(",")
                    AF = ef[0]
                AF_list.append(AF)
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
                genes.append(gene)
                if "frameshift" in consequence or "start_lost" in consequence or "stop_gained" in consequence or "stop_lost" in consequence:
                    lof.append("y")
                else:
                    lof.append("n")
        for i in range(len(genes)):
            print(genes[i], lof[i], file = f)

    count = 0
    for i in range(len(n_shared)):
        count += int(n_shared[i])
    print("Mean number of breeds each variant is shared by: ", count/len(n_shared))
    print("Max number of breeds: ", max(n_shared))
    print("Min number of breeds: ", min(n_shared))

    AF_count = 0
    for i in range(len(AF_list)):
        AF_count += float(AF_list[i])
    print("Mean allele frequency: ", AF_count/len(AF_list))

    print("number of genes affected: ", len(set(genes)))
    print("Number of lof variants: ", lof.count("y"))

    #Get breeds that tend to have the most discrepancies.
    breeds = []
    for key, value in cp_dict.items():
        a = value.split(",")
        for i in range(len(a)):
            breeds.append(a[i])
    print(len(breeds))

    with open("breed_differences.txt", "w") as f:
        breeds_set = list(set(breeds))
        print("Breed_rare_common\tn_variants", file = f)
        for item in range(len(breeds_set)):
            print(breeds_set[item], breeds.count(breeds_set[item]), sep = "\t", file =f)

