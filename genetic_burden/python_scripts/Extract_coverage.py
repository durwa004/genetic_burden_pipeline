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
            help="Path to dir containing the DOC files [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    DOC = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/metrics/DOC_totals/DOC_by_horse.txt"

    with open(DOC, "w") as DOC_file:
        count =0
        overall_DOC = 0
        print("HORSE\ttotal_DOC\tnuclear_placed_DOC", file = DOC_file)
        for file_name in os.listdir(data):
            count +=1
            a = file_name.split("_coverage.tsv")
            horse = a[0]
            print(f"Processing {horse}")
            with open(data + "/" + file_name,"r") as f:
                horse_DOC = 0
                chrom = 0
                for line in f:
                    line = line.rstrip("\n").split("\t")
                    if line[0] == "Total Average:":
                        total = line[1]
                    elif line[0] == "NC_001640.1":
                        continue
                    elif "NC_" in line[0]:
                        chrom +=1
                        b = float(line[1])
                        horse_DOC = horse_DOC + b
                av_DOC = horse_DOC/chrom
                print(horse,total,av_DOC, file = DOC_file, sep = "\t")
            overall_DOC = overall_DOC + av_DOC
            DOC_DOC = overall_DOC/count
            print(DOC_DOC)
