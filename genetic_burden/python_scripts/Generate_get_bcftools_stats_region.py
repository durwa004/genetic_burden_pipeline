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
            help="vcf to split to get variant details from [required]")
    parser.add_argument(
            "-l", "--list_c",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="List of chromosomes to extract information from in .bed format (chrom/start/end) [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    input_file = args.vcf
    chromosomes = args.list_c

#One job for each chromosome (and delete intermediate files - except for top 10 and bottom 10 for each chromosome)

    chrom = {}
    with open(chromosomes) as v_list:
        for line in v_list:
            line = line.rstrip("\n").split("\t")
            chrom[line[0]] = line[2]

    for chromosome in chrom.keys():
        header = (
                  "#!/bin/bash -l\n"  
                  "#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=4g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  f"#PBS -o $PBS_JOBID.Extract_variants_{chromosome}.out\n"
                  f"#PBS -e $PBS_JOBID.Extract_variants_{chromosome}.err\n"
                  f"#PBS -N Extract_variants_{chromosome}.pbs\n"
                 "#PBS -q small\n"
                 "module load bcftools\n\n"
                 )

        pbs = f"Extract_variants_{chromosome}.pbs"
        regions = []
        number = 0
        for i in range(int(chrom[chromosome])):
            if number < (int(chrom[chromosome])-10001):
                a = str(number) + ":" + str((number + 10000))
                regions.append(a)
                number += 10001
            elif number < (int(chrom[chromosome])):
                a = str(number) + ":" + str(chrom[chromosome])
                regions.append(a)
                number += 10001
            else:
                break
 
        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            for i,v in enumerate(regions):
                b = v.split(":")
                print(f"bcftools view -r {chromosome}:", b[0], "-", b[1], f" {input_file} | bcftools stats > {chromosome}_", b[0], "_", b[1], ".stats", sep = "", file = f)
