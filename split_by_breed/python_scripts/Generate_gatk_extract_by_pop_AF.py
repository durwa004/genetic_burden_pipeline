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
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    breeds = ["QH", "TB", "STB", "Shetland", "Morgan","Icelandic", "Clydesdale", "Belgian", "Arabian", "WP"]
    for i in range(len(breeds)):
        breed = breeds[i]
        header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=8,walltime=12:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.gatk_selectvariants_{breed}.out\n"
              f"#PBS -e $PBS_JOBID.gatk_selectvariants_{breed}.err\n"
              f"#PBS -N bcftools_view_{breed}.pbs\n"
              "#PBS -q batch\n"
              "source /home/mccuem/shared/.local/conda/bin/activate gatk4_4.1.0\n"
             )
    
        pbs = os.path.join(os.getcwd(),f"gatk_selectvariants_{breed}.pbs")
    
        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print(f'gatk SelectVariants -R ../../../GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fna -select "AF.0 > 0.10" -V thesis_intersect_{breed}_rare.vcf.gz -O rare_{breed}_common_pop.vcf.gz', file=f)
            print(f'#gatk SelectVariants -R ../../GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fna -select "AF.0 <0.005" -V thesis_intersect_{breed}_common.vcf.gz -O common_{breed}_rare_pop.vcf.gz', file = f)
