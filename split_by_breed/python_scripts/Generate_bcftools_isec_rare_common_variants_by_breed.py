#!/usr/bin/env python3
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
    parser.add_argument(
            "-b", "--breed",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Breed group to extract [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    breed = args.breed

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=8,walltime=04:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.bcftools_isex_{breed}.out\n"
              f"#PBS -e $PBS_JOBID.bcftools_isex_{breed}.err\n"
              f"#PBS -N bcftools_isec_{breed}.pbs\n"
              "#PBS -q batch\n"
              "module load bcftools\n"
             )
    
    pbs = os.path.join(os.getcwd(),f"bcftools_isec_{breed}.pbs")
    breed_list = "Arabian", "Belgian", "Clydesdale", "Icelandic", "Morgan", "QH", "STB", "TB", "WP", "Shetland"

    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        for i in range(len(breed_list)):
            if breed == breed_list[i]:
                continue
            print(f"bcftools isec -p {breed}_rare_",breed_list[i],f"_common {breed}_rare.vcf.gz ",breed_list[i],"_common.vcf.gz", file = f, sep = "")
