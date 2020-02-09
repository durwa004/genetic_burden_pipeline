#!/usr/bin/env python3
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
            help="Path to dir for output files  [required]")
    parser.add_argument(
            "-v", "--vcf",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="annovar output file to extract variants from [required]")
    parser.add_argument(
            "-l", "--list_v",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="List of variants (no header) chrom \t pos format [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    input_file = args.vcf
    variants = args.list_v

    chrom_pos = {}
    with open(variants) as v_list:
        for line in v_list:
            line = line.rstrip("\n").split("\t")
            ab = line[0] + ":" + line[1]
            chrom_pos[ab] = line[2]

#    with open(input_file) as in_file, open(data + "/lof_annovar.txt", "w") as output_file:
    with open(input_file) as in_file, open(data + "/genetic_burden_annovar.txt", "w") as output_file:
        a = 0
        for line in in_file:
            line = line.rstrip("\n").split("\t")
            if len(line) >1:
                b = line[0] + ":" + line[1]
                if a == b:
                    next
                else:
                    a = line[0] + ":" + line[1]
                    if a in chrom_pos.keys():
                       # if line[15] == "y":
                            print("\t".join(line), file = output_file)
                       # else:
                       #     if chrom_pos[a] == "se_only": 
                       #         print("\t".join(line), file = output_file) 
