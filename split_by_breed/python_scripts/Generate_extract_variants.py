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
            help="Path to dir for input chrom/pos files  [required]")
    parser.add_argument(
            "-v", "--vcf",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="vcf to extract variants from [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    input_file = args.vcf

    for filename in os.listdir(data):
        b = filename.split(".txt")
        a = b[0]
        with open(f"Extract_variants_"+ a+".pbs", "w") as f:
                print("#!/bin/bash -l\n"
                  "#PBS -l nodes=1:ppn=8,walltime=06:00:00,mem=1g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  f"#PBS -o $PBS_JOBID.Extract_variants_"+ a + ".out\n"
                  f"#PBS -e $PBS_JOBID.Extract_variants_"+ a+".err\n"
                  f"#PBS -N Extract_variants_"+ a+".pbs\n"
                  "#PBS -q batch\n"
                  "module load bcftools\n"
                  f"cd {data}/../breed_rare_common_snpeff\n"
                  f"bcftools view -R {data}/{filename} {input_file} > {a}_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {a}_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {a}_snpeff.vcf.gz", file = f)

