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
            help="Path to dir containing the annovar output_files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
   
    tmp = []
    for file_name in os.listdir(data):
        if file_name.endswith("_annovar_intersect.variant_function"):
            a = file_name.split("_annovar_intersect.variant_function")
            tmp.append(a[0])

    for i in range(len(tmp)):
        pbs = os.path.join(os.getcwd(), "ANNOVAR_variant_type_" + tmp[i] + ".pbs")
        header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=8,walltime=24:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.ANNOVAR_variant_type_" + tmp[i] + ".out\n"
              "#PBS -e $PBS_JOBID.ANNOVAR_variant_type_" + tmp[i] + ".err\n"
              "#PBS -N ANNOVAR_variant_type_" + tmp[i] + ".pbs\n"
              "#PBS -q small\n"
              "source /home/mccuem/shared/.local/conda/bin/activate snakemake\n"
                 )

        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print(f"python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/extract_annovar_variants_all.py -c ", tmp[i], file = f, sep = "")
