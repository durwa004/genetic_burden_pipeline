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
            help="Path to dir containing the ibio output files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    tmp = []
    for file_name in os.listdir(data + "/../"):
        if file_name.endswith(".genotyped.vcf.gz"):
            a = file_name.split(".genotyped.vcf.gz")
            tmp.append(a[0])

    for i in tmp:
        pbs = os.path.join(os.getcwd(), f"extract_prze_" + i + ".pbs")
        header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=8,walltime=12:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.extract_prze_" + i + ".out\n"
              "#PBS -e $PBS_JOBID.extract_prze_" + i + ".err\n"
              "#PBS -N extract_prze_" + i + ".pbs\n"
              "#PBS -q batch\n"
                 )
        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print("module load bcftools", file = f)
            print("bcftools view -S ../../horse_genomes_IDs.txt --min-ac 1 ../" + i + ".genotyped.vcf.gz > " + i + ".genotyped.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip " + i + ".genotyped.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix " + i + ".genotyped.vcf.gz", file=f)
