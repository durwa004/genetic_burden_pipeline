#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:04:17 2019

@author: jonahcullen
"""

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
            help="Path to dir containing the bcftools stats output files [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    AFR = 0
    AFR_n = 0
    AFC = 0
    AFC_n = 0
    with open(data + "/breed_pop_rare_common_number_of_variants_unique.txt", "w") as info_file, open("bcftools_query_pop_chrom_pos.sh", "w") as output_f, open(data + "/breed_rare_pop_common_number_of_variants_shared", "w") as info_file2, open(data + "/breed_common_pop_rare_number_of_variants_shared", "w") as info_file3:
#    with open(data + "/breed_pop_rare_common_number_of_variants_shared.txt", "w") as info_file, open("bcftools_query_pop_chrom_pos_shared.sh", "w") as output_f, open(data + "/breed_pop_number_of_variants_shared", "w") as info_file2:
        print("Breed\tno_samples\tno_records\tno_SNPs\tno_MNPs\tno_indels\tno_others\tno_multiallelic_sites\tno_nultiallelic_SNPs\tts\ttv\ttstv", file= info_file)
        print("Breed\tno_samples\tno_records\tno_SNPs\tno_MNPs\tno_indels\tno_others\tno_multiallelic_sites\tno_nultiallelic_SNPs\tts\ttv\ttstv", file= info_file2)
        print("Breed\tno_samples\tno_records\tno_SNPs\tno_MNPs\tno_indels\tno_others\tno_multiallelic_sites\tno_nultiallelic_SNPs\tts\ttv\ttstv", file= info_file3)
        for entry in os.listdir(data):
            if os.path.isdir(data + "/" + entry) == True:
                for file_name in os.listdir(data + "/" +entry):
                    if "unique_variants" in entry:
                        if "0000.vcf.stats" in file_name:
                            with open(data + "/" + entry + "/" + file_name,"r") as f:
                                for line in f:
                                    line = line.rstrip("\n").split("\t")
                                    if "#" in line[0]:
                                        next
                                    else:
                                        if line[0] == "ID":
                                            a = line[2].split(".")
                                            chrom = a[0]
                                        elif line[2] == "number of samples:":
                                            sn = line[3]
                                        elif line[2] == "number of records:":
                                            rn = line[3]
                                        elif line[2] == "number of SNPs:":
                                            snp = line[3]
                                        elif line[2] == "number of MNPs:":
                                            mnp = line[3]
                                        elif line[2] == "number of indels:":
                                            indel = line[3]
                                        elif line[2] == "number of others:":
                                            other = line[3]
                                        elif line[2] == "number of multiallelic sites:":
                                            ma = line[3]
                                        elif line[2] == "number of multiallelic SNP sites:":
                                            ma_snp = line[3]
                                        elif line[0] == "TSTV":
                                            ts = line[2]
                                            tv = line[3]
                                            tstv = line[4]
                            print(entry, sn, rn, snp, mnp, indel, other, ma, ma_snp, ts, tv, tstv, file = info_file, sep = "\t")
                        elif file_name == "0000.vcf.gz":
                            print(f"cd {data}/{entry}\n", file=output_f)
                            print(f"bcftools query -f '%CHROM\ xxxxxt%POS\ xxxxxn' {file_name} > {data}/breed_pop_rare_common_chrom_pos/{entry}_unique.txt", file = output_f)
                    elif "pop_rare" in entry:
                        if "0002.vcf.stats" in file_name:
                            with open(data + "/" + entry + "/" + file_name,"r") as f:
                                for line in f:
                                    line = line.rstrip("\n").split("\t")
                                    if "#" in line[0]:
                                        next
                                    else:
                                        if line[0] == "ID":
                                            a = line[2].split(".")
                                            chrom = a[0]
                                        elif line[2] == "number of samples:":
                                            sn = line[3]
                                        elif line[2] == "number of records:":
                                            rn = line[3]
                                        elif line[2] == "number of SNPs:":
                                            snp = line[3]
                                        elif line[2] == "number of MNPs:":
                                            mnp = line[3]
                                        elif line[2] == "number of indels:":
                                            indel = line[3]
                                        elif line[2] == "number of others:":
                                            other = line[3]
                                        elif line[2] == "number of multiallelic sites:":
                                            ma = line[3]
                                        elif line[2] == "number of multiallelic SNP sites:":
                                            ma_snp = line[3]
                                        elif line[0] == "TSTV":
                                            ts = line[2]
                                            tv = line[3]
                                            tstv = line[4]
                            print(entry, sn, rn, snp, mnp, indel, other, ma, ma_snp, ts, tv, tstv, file = info_file2, sep = "\t")
                        elif file_name == "0003.vcf.gz":
                            with gzip.open(data + "/" + entry + "/" + file_name,"rt") as f:
                                for line in f:
                                    line = line.rstrip("\n").split("\t")
                                    if "#" in line[0]:
                                        next
                                    else:
                                        info = line[7].split("AF=")
                                        info = info[1].split(";")
                                        info = info[0].split(",")
                                        AFR += float(info[0])
                                        AFR_n +=1
                            print(f"cd {data}/{entry}\n", file=output_f)
                            print(f"bcftools query -f '%CHROM\ xxxxxt%POS\ xxxxxn' {file_name} > {data}/breed_pop_rare_common_chrom_pos/{entry}_shared.txt", file = output_f)
                    elif "pop_common" in entry:
                        if "0003.vcf.stats" in file_name:
                            with open(data + "/" + entry + "/" + file_name,"r") as f:
                                for line in f:
                                    line = line.rstrip("\n").split("\t")
                                    if "#" in line[0]:
                                        next
                                    else:
                                        if line[0] == "ID":
                                            a = line[2].split(".")
                                            chrom = a[0]
                                        elif line[2] == "number of samples:":
                                            sn = line[3]
                                        elif line[2] == "number of records:":
                                            rn = line[3]
                                        elif line[2] == "number of SNPs:":
                                            snp = line[3]
                                        elif line[2] == "number of MNPs:":
                                            mnp = line[3]
                                        elif line[2] == "number of indels:":
                                            indel = line[3]
                                        elif line[2] == "number of others:":
                                            other = line[3]
                                        elif line[2] == "number of multiallelic sites:":
                                            ma = line[3]
                                        elif line[2] == "number of multiallelic SNP sites:":
                                            ma_snp = line[3]
                                        elif line[0] == "TSTV":
                                            ts = line[2]
                                            tv = line[3]
                                            tstv = line[4]
                            print(entry, sn, rn, snp, mnp, indel, other, ma, ma_snp, ts, tv, tstv, file = info_file3, sep = "\t")
                        elif file_name == "0003.vcf.gz":
                            with gzip.open(data + "/" + entry + "/" + file_name,"rt") as f:
                                for line in f:
                                    line = line.rstrip("\n").split("\t")
                                    if "#" in line[0]:
                                        next
                                    else:
                                        info = line[7].split("AF=")
                                        info = info[1].split(";")
                                        info = info[0].split(",") 
                                        AFC += float(info[0])
                                        AFC_n +=1
                            print(f"cd {data}/{entry}\n", file=output_f)
                            print(f"bcftools query -f '%CHROM\ xxxxxt%POS\ xxxxxn' {file_name} > {data}/breed_pop_rare_common_chrom_pos/{entry}_shared.txt", file = output_f)
    print("AF population common = ", AFC/AFC_n)
    print("AF population rare = ", AFR/AFR_n)
