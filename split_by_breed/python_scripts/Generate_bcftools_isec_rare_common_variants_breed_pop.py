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
            help="List of horse IDs and horse breeds [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    horse_ids = args.ids

    breed={}
    with open(horse_ids) as ids:
        ids.readline()
        for line in ids:
            horse,group = line.rstrip("\n").split("\t")
            if group == "Other":
                pass
            else:
                if group in breed.keys():
                    a = breed[group] + ":" + horse
                    breed[group] = a
                else:
                    breed[group] = horse

    for key,value in breed.items():
        with open(f"bcftools_isec_{key}_pop.pbs", "w") as f:
            horses = value.split(":")
            print("#!/bin/bash -l\n"
                "#PBS -l nodes=1:ppn=8,walltime=12:00:00,mem=2g\n"
                "#PBS -m abe\n"
                "#PBS -M durwa004@umn.edu\n"
                f"#PBS -o $PBS_JOBID.bcftools_isec_{key}_pop.out\n"
                f"#PBS -e $PBS_JOBID.bcftools_isec_{key}_pop.err\n"
                f"#PBS -N bcftools_isec_{key}_pop.pbs\n"
                "#PBS -q small\n"
                "module load bcftools\n", file = f)
            print(f"cd {data}\n", file=f)
            #Rare breed/common pop
            print(f"bcftools isec -p {key}_rare_pop_common breed_rare_common_vcfs/{key}_rare.vcf.gz without_{key}_common.vcf.gz", file =f)
            #Common breed/rare pop
            #print(f"bcftools isec -p {key}_common_pop_rare ../breed_rare_common_vcfs/{key}_common.vcf.gz without_{key}_rare.vcf.gz", file =f)
            #Present in breed/absent in pop
            #print(f"bcftools isec -p unique_variants_{key} ../thesis_intersect_{key}.vcf.gz thesis_intersect_without_{key}.vcf.gz", file =f)
