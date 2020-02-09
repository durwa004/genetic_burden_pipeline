#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:04:17 2019
@author: jonahcullen
"""

import argparse
import os
import gzip

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to dir for output files  [required]")
    parser.add_argument(
            "-v", "--vcf",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="vcf to split to get variant details from [required]")
    parser.add_argument(
            "-l", "--list_r",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="List of regions to extract information from [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    input_file = args.vcf
    regions = args.list_r

    filename = regions.split("/")
    filename = filename[-1].split(".txt")
    filename = filename[0]
#One job for each chromosome (and delete intermediate files - except for top 10 and bottom 10 for each chromosome)

    with open(regions) as v_list, open(data + f"/{filename}_regions.txt", "w") as f:
        v_list.readline()
        for line in v_list:
            line = line.rstrip("\n").split("\t")
            print(line[1], line[2], line[3], sep = "\t", file = f)

    header = (
                  "#!/bin/bash -l\n"  
                  "#PBS -l nodes=1:ppn=1,walltime=12:00:00,mem=4g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  f"#PBS -o $PBS_JOBID.Extract_variants_{filename}.out\n"
                  f"#PBS -e $PBS_JOBID.Extract_variants_{filename}.err\n"
                  f"#PBS -N Extract_variants_{filename}.pbs\n"
                 "#PBS -q small\n"
                 "module load bcftools\n\n"
                 )

    pbs = f"Extract_variants_{filename}.pbs"
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {data}\n", file=f)
        print(f"bcftools view -R {filename}_regions.txt {input_file} -o {filename}_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {filename}_snpeff.vcf &&  /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {filename}_snpeff.vcf.gz && java -Djava.io.tmpdir=/home/mccuem/shared/Projects/HorseGenomeProject/Data/temp_files -Xmx12g -jar /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/variant_annotation/SnpEff/snpEff/SnpSift.jar filter ", '"((ANN[0].IMPACT has ',"'HIGH'",') | (ANN[*].IMPACT = ',"'MODERATE'", ') | (ANN[*].IMPACT = ',"'LOW'))",'"', f" {filename}_snpeff.vcf.gz > {filename}_snpeff.coding.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {filename}_snpeff.coding.vcf &&  /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {filename}_snpeff.codingvcf.gz", sep = "", file = f)
