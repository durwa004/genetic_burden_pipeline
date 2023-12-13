#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:04:17 2019

@author: jonahcullen
"""

import argparse
import os
def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to dir containing the bcftools stats output file [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

#Add in breed info:
    breed = {}
    with open("horse_genomes_breeds_tidy.txt", "r") as input_f:
        input_f.readline()
        for line in input_f:
            line = line.rstrip("\n").split("\t")
            breed[line[0]] = line[1]

    with open(data, "r") as input_f, open("ind_number_of_variants.txt", "w") as info_file:
        print("breed\tID\tnRefHom\tnNonRefHom\tnHets\tnSNPs\tnTranstions\tnTransversions\ttstv\tindels\tDOC\tnSingletons\tnMissing", sep = "", file = info_file)
        for line in input_f:
            line = line.rstrip("\n").split("\t")
            if "PSC" == line[0]:
                variants = int(line[4]) + int(line[5])
                tstv = int(line[6])/int(line[7])
                print(breed[line[2]], "\t".join(line[2:6]), variants, "\t".join(line[6:8]), tstv, "\t".join(line[8:11]), line[-1], sep = "\t", file = info_file)
