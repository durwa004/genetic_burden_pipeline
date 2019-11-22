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

    DOC = f"{data}/summary_files/DOC_by_horse.txt"
    reads = f"{data}/summary_files/reads_by_horse.txt"

    with open(DOC, "w") as DOC_file, open(reads, "w") as reads_file:
        count =0
        overall_DOC = 0
        print("HORSE\ttotal_DOC\tnuclear_placed_DOC", file = DOC_file)
        for directory in os.listdir(data):
            for file_name in os.listdir(data + "/" + directory):
                if file_name == "coverage.tsv":
                    count +=1
                    a = directory.split("_")
                    horse = a[0]
                    print(f"Processing {horse}")
                    with open(data + "/" + directory + "/" + file_name,"r") as f:
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
                elif file_name.endswith(".alignment_summary_metrics"):
                    with open(data + "/" + directory + "/" + file_name, "r") as f2:
                        b = file_name.split(".metrics.alignment_summary_metrics")
                        horse = b[0]
                        for line in f2:
                            line = line.rstrip("\n").split("\t")
                            if line[0] == "PAIR":
                                print(horse, line[15], line[16], file = reads_file, sep = "\t")
                        
