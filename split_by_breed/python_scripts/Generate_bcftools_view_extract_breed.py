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
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    horse_ids = os.path.abspath(args.ids)

    breed={}
    with open(horse_ids) as ids:
        ids.readline()
        for line in ids:
            horse,group = line.rstrip("\n").split("\t")
            if group in breed.keys():
                a = breed[group] + ":" + horse
                breed[group] = a
            else:
                breed[group] = horse

    breed_groups = list(set(breed.keys()))

    for i,v in enumerate(breed_groups):
        with open(f"{data}/remove_{v}_ids.list", "w") as output_file:
            for key,value in breed.items():
                if key == v:
                    continue
                else:
                    horses = value.split(":")
                    for i in range(len(horses)):
                        print(horses[i], file = output_file)

    for key,value in breed.items():
        with open(f"bcftools_view_remove_{key}.pbs", "w") as f:
            print("#!/bin/bash -l\n"
                "#PBS -l nodes=1:ppn=8,walltime=12:00:00,mem=2g\n"
                "#PBS -m abe\n"
                "#PBS -M durwa004@umn.edu\n"
                f"#PBS -o $PBS_JOBID.bcftools_view_{key}.out\n"
                f"#PBS -e $PBS_JOBID.bcftools_view_{key}.err\n"
                f"#PBS -N bcftools_view_{key}.pbs\n"
                "#PBS -q batch\n"
                "module load bcftools\n", file = f)
            print(f"cd {data}\n", file=f)
            print(f"#bcftools view -S remove_{key}_ids.list --min-ac 1 ../joint_intersect_without_Prze/thesis_intersect.vcf.gz > thesis_intersect_without_{key}.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_intersect_without_{key}.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_intersect_without_{key}.vcf.gz",file=f) 
