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
            help="Text file containing horse ids and dz - one per line [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    horse_ids = os.path.abspath(args.ids)

    breed = []
    with open(horse_ids) as ids:
        for line in ids:
            horse,group = line.rstrip("\n").split("\t")
            a = group + ":" + horse
            breed.append(a)

    for i,v in enumerate(breed):
        ab = v.split(":")
        ID = ab[0]
        dz = ab[1]
        with open(f"bcftools_view_{ID}_{dz}.pbs", "w") as f:
            print("#!/bin/bash -l\n"
                "#PBS -l nodes=1:ppn=1,walltime=24:00:00,mem=2g\n"
                "#PBS -m abe\n#PBS -M durwa004@umn.edu\n"
                f"#PBS -o $PBS_JOBID.bcftools_view_{ID}_{dz}.out\n"
                f"#PBS -e $PBS_JOBID.bcftools_view_{ID}_{dz}.err\n"
                f"#PBS -N bcftools_view_{ID}_{dz}.pbs\n"
                "#PBS -q small\n"
                "module load bcftools\n", file = f)
            print(f"cd {data}\n", file=f)
            print(f"#bcftools view -s {dz} --min-ac 1 ../../../joint_union/thesis_union.vcf.gz > thesis_union_{ID}_{dz}.vcf && \n/home/mccuem/durwa004/bin/vt/vt decompose -s thesis_union_{ID}_{dz}.vcf -o thesis_union_{ID}_{dz}.decomposed.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_union_{ID}_{dz}.decomposed.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_union_{ID}_{dz}.decomposed.vcf.gz && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_union_{ID}_{dz}.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_union_{ID}_{dz}.vcf.gz",file=f)    
