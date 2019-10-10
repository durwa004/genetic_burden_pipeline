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
            "-i", "--ids",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Text file containing horse ids and breed - one per line [required]")
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
    horse_ids = os.path.abspath(args.ids)
    breed = args.breed

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=12,walltime=48:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.bcftools_view_{breed}.out\n"
              f"#PBS -e $PBS_JOBID.bcftools_view_{breed}.err\n"
              f"#PBS -N bcftools_view_{breed}.pbs\n"
              "#PBS -q small\n"
             )
    
    out = f"{data}/{breed}_ids.list")
    with open(horse_ids) as ids, open(out, "w") as output_file:
        for line in ids:
            horse,group = line.rstrip("\n").split("\t")
            if group == breed:
                print(horse, file = output_file)
    
    pbs = os.path.join(os.getcwd(),f"bcftools_view_{breed}.pbs")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        count = 0
        print(f"bcftools stats -s {breed}_ids.list thesis_intersect.vcf.gz > thesis_intersect_{breed}.stats",file=f) 
        print(f"Creating .pbs script for {breed}\nNumber of horses: {count}")
