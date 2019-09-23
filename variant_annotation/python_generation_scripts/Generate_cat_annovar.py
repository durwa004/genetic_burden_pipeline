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
            help="Path to dir containing the annovar output files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=1,walltime=08:00:00,mem=2g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.annovar_cat.out\n"
              f"#PBS -e $PBS_JOBID.annovar_cat.err\n"
              f"#PBS -N annovar_cat.pbs\n"
              "#PBS -q small\n"
             )
    
    tmp = []
    for file_name in os.listdir(data):
        if "NC_" in file_name:
            if file_name.endswith("_annovar_intersect.exonic_variant_function"):
                tmp.append(file_name)
    
    pbs = os.path.join(os.getcwd(), "annovar_cat.pbs")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        print("cat", " ".join(tmp), " > annovar_coding.txt",file=f)
