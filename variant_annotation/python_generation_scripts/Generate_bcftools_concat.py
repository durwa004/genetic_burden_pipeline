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
            help="Path to dir containing the annovar/SnpEff union output files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=12,walltime=48:00:00,mem=12g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.snpeff_concat.out\n"
              f"#PBS -e $PBS_JOBID.snpeff_concat.err\n"
              f"#PBS -N snpeff_concat.pbs\n"
              "#PBS -q small\n"
              "module load bcftools\n"
             )
    
    tmp = []
    for file_name in os.listdir(data):
        if file_name.endswith(".coding.ann.vcf.gz"):
            tmp.append(file_name)
    
    pbs = os.path.join(os.getcwd(), "snpeff_concat.pbs")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        print("bcftools concat", " ".join(tmp), " -o SnpEff_coding.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip SnpEff_coding.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix SnpEff_coding.vcf.gz",file=f)
