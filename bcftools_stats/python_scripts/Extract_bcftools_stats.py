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
    parser.add_argument(
            "-p", "--prog",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Program that we are getting the bcfstats info from (i.e. bcftools/gatk/union/intersect [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    program = args.prog

    info = data + f"/{program}_number_of_variants.txt"

    with open(info, "w") as info_file:
        print("CHROM\tno_samples\tno_records\tno_SNPs\tno_MNPs\tno_indels\tno_others\tno_multiallelic_sites\tno_nultiallelic_SNPs\tts\ttv\ttstv", file= info_file)
        for file_name in os.listdir(data):
            if file_name.endswith(".stats"):
                print(f"Processing {file_name}")
                with open(data + "/" + file_name,"r") as f:
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
                    print(chrom, sn, rn, snp, mnp, indel, other, ma, ma_snp, ts, tv, tstv, file = info_file, sep = "\t")
