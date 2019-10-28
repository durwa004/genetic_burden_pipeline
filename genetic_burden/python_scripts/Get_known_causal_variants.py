known = {}
with open("known_SNVs.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        known[line[3]] = line[2]

known_checked = {}
with open("known_dz_variants.txt", "r") as input_file, open("known_variants_AFs.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] in known.keys():
            if known[line[1]] == line[0]:
                known_checked[line[1]] = line[0]
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                AC = bc[1]
                cd = ab[1].split("AF=")
                AF = cd[1]
                if "," in AF:
                    ef = AF.split(",")
                    AF = ef[0]
                else:
                    next
                print("\t".join(line[:5]), AC, AF, "\t".join(line[9:]), file = output_file, sep = "\t")

with open("known_dz_variants.txt", "r") as input_file, open("known_variants_AFs.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if line[1] in known_checked.keys():
            if known_checked[line[1]] == line[0]:
                ab = line[7].split(";")
                bc = ab[0].split("AC=")
                AC = bc[1]
                cd = ab[1].split("AF=")
                AF = cd[1]
                if "," in AF:
                    ef = AF.split(",")
                    AF = ef[0]
                else:
                    AF = cd[1]
                print("\t".join(line[:5]), AC, AF, "\t".join(line[9:]), file = output_file, sep = "\t")

#Get list of horse ids in order of vcf.
import gzip
with gzip.open("../thesis_intersect.vcf.gz", "rt") as input_file, open("horse_ids_vcf.txt", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#CHROM" in line[0]:
            for i in range(len(line)):
                print(line[i], file = output_file)


