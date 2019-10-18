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
    parser.add_argument(
            "-m3", "--mx3",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="3% AF AC cut off (unique for each breed) [required]")
    parser.add_argument(
            "-m5", "--mn5",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="5% AF AC cut off (unique for each breed) [required]")            
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    breed = args.breed
    mx = args.mx3
    mn = args.mn5

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=8,walltime=12:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.bcftools_view_{breed}_mx_AC_{mx}.out\n"
              f"#PBS -e $PBS_JOBID.bcftools_view_{breed}_mx_AC_{mx}.err\n"
              f"#PBS -N bcftools_view_{breed}_mx_AC_{mx}.pbs\n"
              "#PBS -q batch\n"
              "module load bcftools\n"
             )
    
    pbs = os.path.join(os.getcwd(),f"bcftools_view_{breed}_mx_AC_{mx}.pbs")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        print(f"bcftools view --max-ac {mx} --min-af 0.1 thesis_intersect_{breed}.vcf.gz > rare_{breed}_common_population.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip rare_{breed}_common_population.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix rare_{breed}_common_population.vcf.gz", file = f, sep = "")
        print(f"#bcftools view --min-ac {mn} --max-af 0.005 thesis_intersect_{breed}.vcf.gz > common_{breed}_rare_population.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip common_{breed}_rare_population.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix common_{breed}_rare_population.vcf.gz", file = f, sep = "")
