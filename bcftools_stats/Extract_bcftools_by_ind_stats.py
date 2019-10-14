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
                print(f"Processing {file_name}")
                AF = []
                freq = []
                with open(data + "/" + file_name,"r") as f:
                    for line in f:
                        line = line.rstrip("\n").split("\t")
                        if "#" in line[0]:
                            next
                        else:
                             if line[0] == "PSC":
                                 sample = line[2]
                                 nRefHom = line[3]
                                 nNonRefHom = line[4]
                                 nHets = line[5]
                                 nTs = line[6]
                                 nTv = line[7]
                                 nIndels = line[8]
                                 nSingletons = line[10]
                             elif line[0] == "HWE":
                                 AF.append(line[2])
                                 freq.append(line[3])
                print(sample,breed,nRefHom,nNonRefHom,nHets,nTs,nTv,nIndels,nSingletons, file = info_file, sep = "\t")
                with open(f"{data}/{sample}_{breed}_AF_freq.txt", "w") as output_file:
                    print("AF\tfrequency",file=output_file)
                    for item in range(len(AF)):
                        print(AF[item], freq[item], file = output_file, sep = "\t") 
