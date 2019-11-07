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

    info = data + "/intersect_by_ind_number_of_variants.txt"

    with open(info, "w") as info_file:
        print("Sample\tbreed\tnRefHom\tnNonRefHom\tnHets\tTs\tTv\tnIndels\tnSingletons", file= info_file)
        for file_name in os.listdir(data):
            if file_name.endswith(".stats"):
                a = file_name.split("_")
                b = a[1].split(".stats")
                breed = b[0]
                sample = a[0]
                print(f"Processing {file_name}")
                AF = []
                freq = []
                with open(data + "/" + file_name,"r") as f, open(f"{data}/{sample}_{breed}_AF_freq.txt", "w") as output_file:
                    print("AF\tfrequency",file=output_file)
                    for line in f:
                        line = line.rstrip("\n").split("\t")
                        if "#" in line[0]:
                            next
                        elif line[0] == "PSC":
                             print(line[2],breed,line[3],line[4],line[5],line[6],line[7],line[8],line[10],file=info_file,sep="\t")
                        elif line[0] == "HWE":
                             AF.append(line[2])
                             freq.append(line[3])
                    for item in range(len(AF)):
                        print(AF[item], freq[item], file = output_file, sep = "\t") 
