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
            help="Path to dir containing the ibio intersect/union by chromosome files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
   
    temp = "-Djava.io.tmpdir=/home/mccuem/shared/Projects/HorseGenomeProject/Data/temp_files"
    snpeff = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/variant_annotation/SnpEff/snpEff/"
    version = "EquCab3"
 
    print(f"Using SnpEff version {version}")

    tmp = []
    for file_name in os.listdir(data):
        if file_name.endswith(".vcf.gz"):
            a = file_name.split("_intersect_joint.vcf.gz")
            tmp.append(a[0])

    for i in tmp:
        pbs = os.path.join(os.getcwd(), "SnpEff_intersect_" + i + ".pbs")
        header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=4,walltime=24:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.SnpEff_intersect_" + i + ".out\n"
              "#PBS -e $PBS_JOBID.SnpEff_intersect_" + i + ".err\n"
              "#PBS -N SnpEff_intersect_" + i + ".pbs\n"
              "#PBS -q small\n"
                 )

        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print(f"java {temp} -Xmx4g -jar {snpeff}snpEff.jar -C {snpeff}snpEff.config {version} " + i + "_intersect_joint.vcf.gz > SnpEff/" + i + "_snpeff_intersect.ann.vcf -csvstats SnpEff/" + i + "_snpeff_intersect.csv -stats SnpEff/" + i + "_snpeff_intersect.html", file=f, sep="")
            print("/home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip SnpEff/" + i + "_snpeff_intersect.ann.vcf",file=f,sep="")
            print("/home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix SnpEff/" + i + "_snpeff_intersect.ann.vcf.gz", file=f,sep="") 
