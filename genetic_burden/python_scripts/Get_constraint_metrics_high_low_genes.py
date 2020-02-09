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
            help="Path to dir for output files  [required]")
    parser.add_argument(
            "-t", "--txt",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="vcf to extract variants from [required]")
    return parser

if __name__ == '__main__':

    parser = make_arg_parser()
    args = parser.parse_args()

    path = os.path.abspath(args.data)
    input_vcf = args.txt

#Add in gene constraint information
    HI = {}
    with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/disease_case_analysis/constraint_files/HI_Predictions_Version3.bed", "r") as input_file:
        input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            a = line[3].split("|")
            HI[a[0]] = line[4]

    pLI = {}
    with open("/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/disease_case_analysis/constraint_files/gnomad.v2.1.1.lof_metrics.by_transcript.txt", "r") as input_file:
        a = input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            if line[2] == "true":
                pLI[line[0]] = line[21]

#Goal is to print a list of HI and pLI scores for the genes in the high variation regions and low variation regions
filename = input_vcf.split("/")
filename = filename[-1].split(".txt")
genes = []
with open(input_vcf, "r") as input_file, open(path + "/" + filename[0] + "_constraint.txt", "w") as output_file:
    input_file.readline()
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        genes.append(line[9])
    genes = list(set(genes))
    for i in range(len(genes)):
        if genes[i] in HI.keys():
            HI_i = HI[genes[i]]
        else:
            HI_i = "NA"
        if genes[i] in pLI.keys():
            pLI_i = pLI[genes[i]]
        else:
            pLI_i = "NA"
        print(genes[i], HI_i, pLI_i, file = output_file, sep = "\t")
