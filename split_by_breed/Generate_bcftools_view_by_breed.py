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
    parser.add_argument(
            "-m", "--maxAC",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Maximum AC for rare variants [required]")
    parser.add_argument(
            "-c", "--com",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Minimum AC for common variants [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    horse_ids = os.path.abspath(args.ids)
    breed = args.breed
    maxAC = args.maxAC
    common = args.com

    header = (
              "#!/bin/bash -l\n"  
              "#PBS -l nodes=1:ppn=12,walltime=48:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              f"#PBS -o $PBS_JOBID.bcftools_view_{breed}.out\n"
              f"#PBS -e $PBS_JOBID.bcftools_view_{breed}.err\n"
              f"#PBS -N bcftools_view_{breed}.pbs\n"
              "#PBS -q small\n"
              "module load bcftools\n"
             )
    
    out = f"{data}/{breed}_ids.list"
    count = 0
    with open(horse_ids) as ids, open(out, "w") as output_file:
        for line in ids:
            horse,group = line.rstrip("\n").split("\t")
            if group == breed:
                count +=1
                print(horse, file = output_file)
    
    pbs = os.path.join(os.getcwd(),f"bcftools_view_{breed}.pbs")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        print(f"bcftools view -S {breed}_ids.list --min-ac 1 thesis_intersect.vcf.gz > thesis_intersect_{breed}.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_intersect_{breed}.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_intersect_{breed}.vcf.gz",file=f) 
        print(f"bcftools view thesis_intersect_{breed}.vcf.gz --min-ac 1 --max-ac {maxAC} > thesis_intersect_{breed}_rare.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_intersect_{breed}_rare.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_intersect_{breed}_rare.vcf.gz", file=f)
        print(f"bcftools view thesis_intersect_{breed}.vcf.gz --min-ac {common} > thesis_intersect_{breed}_common.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_intersect_{breed}_common.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_intersect_{breed}_common.vcf.gz", file=f)
        print(f"Creating .pbs script for {breed}\nNumber of horses: {count}")
