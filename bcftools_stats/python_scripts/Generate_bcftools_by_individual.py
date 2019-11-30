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
            help="Path to dir containing the ibio output files of interest [required]")
    parser.add_argument(
            "-ind", "--individual",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="List of individuals to get bcftools to run on [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    horses = os.path.abspath(args.individual) 

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=4,walltime=06:00:00,mem=2g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.bcftools_stats_by_ind.out\n"
              f"#PBS -e $PBS_JOBID.bcftools_stats_by_ind.err\n"
              f"#PBS -N bcftools_stats_by_ind.pbs\n"
              "#PBS -q small\n"
              "module load bcftools\n"
             )
    
    with open(horses) as input_file:
        for line in input_file:
            horse,breed = line.rstrip("\n").split("\t")
            pbs = os.path.join(os.getcwd(), f"bcftools_stats_{horse}.pbs")
            #create job script for each individual horse
            with open(pbs, "w") as f:
                print(header, file=f)
                print(f"cd {data}\n", file=f)
                print(f"bcftools stats -s {horse} thesis_union.vcf.gz > ind_bcftools_stats_files/{horse}_{breed}.stats", file=f) 
