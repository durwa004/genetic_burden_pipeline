#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:04:17 2019
@author: jonahcullen
"""

import argparse
import os
import gzip
from statistics import mean

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to dir for input chrom/pos files  [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    chrom_pos = {}
    for filename in os.listdir(data):
#        if filename.endswith(".txt"):
        if filename.endswith("_shared.txt"):
            b = filename.split(".txt")
            a = b[0]
            with open(data + "/" + filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if len(line) <2:
                        continue
                    else:
                        a = line[0] + ":" + line[1]
                        c = filename.split(".txt")
                        c = c[0]
                    if a in chrom_pos.keys():
                        b = chrom_pos[a]
                        b = b + ":" + c
                        chrom_pos[a] = b
                    else:
                        chrom_pos[a] = c

        #Get average number of times each variant appears
    count = 0
    min_c = 3
    max_c = 0
    for key, value in chrom_pos.items():
        a = value.split(":")
        count += len(a)
        if len(a) < int(min_c):
            min_c = len(a)
        elif len(a) > int(max_c):
            max_c = len(a)
        
    print("Variant details. Mean = ")
    print(count/(len(chrom_pos.keys())))
    print("Max = ")
    print(max_c)
    print("Min = ")
    print(min_c)

#    with open(f"{data}/rare_common_breed_chrom_pos.txt", "w") as f, open(f"{data}/rare_common_breed_chrom_pos_with_breed_details.txt", "w") as f2:
#    with open(f"{data}/rare_common_breed_pop_chrom_pos.txt", "w") as f, open(f"{data}/rare_common_breed_pop_chrom_pos_with_breed_details.txt", "w") as f2:
    with open(f"{data}/shared_breed_pop_chrom_pos.txt", "w") as f, open(f"{data}/shared_pop_chrom_pos_with_breed_details.txt", "w") as f2:
        for key,value in chrom_pos.items():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)
            print(a[0], a[1], value, sep = "\t", file = f2)

    with open(f"Extract_variants_shared_pop.pbs", "w") as f:
        print("#!/bin/bash -l\n"
                  "#PBS -l nodes=1:ppn=1,walltime=12:00:00,mem=1g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  "#PBS -o $PBS_JOBID.Extract_variants_shared_pop.out\n"
                  "#PBS -e $PBS_JOBID.Extract_variants_shared_pop.err\n"
                  "#PBS -N Extract_variants_shared_pop.pbs\n"
                  "#PBS -q small\n"
                  "module load bcftools\n"
                  f"cd {data}/../breed_pop_rare_common_snpeff\n"
                  f"bcftools view -R {data}/shared_pop_chrom_pos.txt ../../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz > shared_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip shared_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix shared_pop_snpeff.vcf.gz", file = f)

