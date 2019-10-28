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
            help="vcf to extract variants from [required]")
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

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=2,walltime=12:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.Extract_variants.out\n"
              "#PBS -e $PBS_JOBID.Extract_variants.err\n"
              "#PBS -N Extract_variants.pbs\n"
              "#PBS -q small\n"
             )

    chrom_pos = []
    impact = []
    with open(variants) as v_list:
        for line in v_list:
            line = line.rstrip("\n").split("\t")
            ab = line[0] + ":" + line[1] + "-" + line[1]
            chrom_pos.append(ab)
            bc = line[2] + "_" + line[3] + "_"
            impact.append(bc)

    pbs = os.path.join(os.getcwd(),"Extract_variants.pbs")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        for i in range(len(chrom_pos)):
            print(f"/home/mccuem/shared/.local/bin/tabix {input_file} ", chrom_pos[i], " > ", impact[i], chrom_pos[i], ".txt", file = f, sep = "")
