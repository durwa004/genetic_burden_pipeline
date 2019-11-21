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
            "-id", "--idlist",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Text file including disease and horse id [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    ids = os.path.abspath(args.idlist)

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
    
    dz = {}
    with open(ids) as input_file:
        input_file.readline()
        for line in input_file:
            line = line.rstrip("\n").split("\t")
            if line[1] in dz.keys():
                dz[line[1]][line[0]] = 0
            else:
                dz[line[1]]= {}
                dz[line[1]][line[0]] = 0
    dz["hypertrophy"]["M10643"] = 0    
    
    for disease in dz.keys():
        pbs = os.path.join(os.getcwd(),f"bcftools_stats_{disease}.pbs")
        with open(pbs, "w") as f:
            print(header, file = f)
            print(f"cd {data}\n", file = f)
            horses = list(dz[disease].keys())
            print(f"bcftools stats -s ", ",".join(horses), f" thesis_intersect_snpeff.ann.vcf.gz > {disease}.stats", file = f, sep = "")
