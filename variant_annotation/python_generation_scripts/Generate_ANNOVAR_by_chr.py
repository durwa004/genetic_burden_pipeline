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
            help="Path to dir containing the ibio intersect/union by chromosome files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
   
    annovar_v = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/variant_annotation/annovar/"
    annovar_db = "GCF_002863925.1_EquCab3.0_db"
 
    print(f"Using {annovar_db} for ANNOVAR version {annovar_v}")

    tmp = []
    for file_name in os.listdir(data):
        if file_name.endswith(".vcf.gz"):
            a = file_name.split("_intersect_joint.vcf.gz")
            tmp.append(a[0])

    for i in tmp:
        pbs = os.path.join(os.getcwd(), "ANNOVAR_intersect_" + i + ".pbs")
        header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=8,walltime=24:00:00,mem=20g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.ANNOVAR_intersect_" + i + ".out\n"
              "#PBS -e $PBS_JOBID.ANNOVAR_intersect_" + i + ".err\n"
              "#PBS -N ANNOVAR_intersect_" + i + ".pbs\n"
              "#PBS -q small\n"
                 )

        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print(f"#perl {annovar_v}convert2annovar.pl -format vcf4 -allsample -include -withfreq " + i + "_intersect_joint.vcf.gz --out annovar/" + i + "_annovar_intersect.avinput\n"
                  + f"perl {annovar_v}annotate_variation.pl -buildver GCF_002863925.1_EquCab3.0 annovar/" + i + f"_annovar_intersect.avinput {annovar_v}{annovar_db} --geneanno --hgvs --aamatrixfile {annovar_v}example/grantham.matrix --out annovar/" + i + "_annovar_intersect", file = f, sep = "")
