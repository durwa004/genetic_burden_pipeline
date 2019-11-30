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
            "-p", "--program",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Program the chromosomes are from [required]")
    parser.add_argument(
            "-i", "--ids",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="List of horse ids to check [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    prog = args.program
    horse_ids = os.path.abspath(args.ids)

    info = data + f"/{prog}_ind_number_of_variants.txt"

    #Create dictionary for horses
    chrom_dict = {}
    with open(horse_ids, "r") as input_file:
        input_file.readline()
        for line in input_file:
             line = line.rstrip("\n").split("\t")
             a = line[0] + "." + line[1]
             chrom_dict[a] = {}

    #Add in chromosome for each horse
    for filename in os.listdir(data + "/../"):
        if filename.endswith(".genotyped.vcf.gz"):
            a = filename.split(".genotyped")
            for horse in chrom_dict.keys():
                b = "0:0:0:0:0:0"
                chrom_dict[horse][a[0]] = b

    for filename in os.listdir(data):
        if filename.endswith(".stats"):
            a = filename.split("_")
            horse = a[0] + "." + a[1]
            if "unplaced" in a[2]:
                c = a[2].split(".stats")
                b = c[0]
            elif a[3] == "001640":
                b = "NC_" + a[3] + "_1"
            else:
                b = "NC_" + a[3] + "_3"
            with open(f"{data}/{filename}") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if line[0] == "PSC":
                        details = line[3] + ":" + line[4] + ":" + line[5] + ":" + line[6] + ":" + line[7] + ":" + line[8]
                        chrom_dict[horse][b] = details 

    with open(info, "w") as info_file:
        print("sample.breed\tchromosome\tnRefHom\tnNonRefHom\tnHets\tTs\tTv\tnIndels", file= info_file)
        for horse in chrom_dict.keys():
            for chrom in chrom_dict[horse]:
                a = chrom_dict[horse][chrom].split(":")
                print(horse,chrom,"\t".join(a), sep = "\t",file = info_file)
