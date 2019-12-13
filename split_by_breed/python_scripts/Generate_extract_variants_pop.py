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
            help="Path to dir for input chrom/pos files  [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    chrom_pos = {}

    for filename in os.listdir(data):
        b = filename.split(".txt")
        a = b[0]
        if "_rare_pop_common_" in b[0]:
            with open(data + "/" + filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if len(line) <2:
                        continue
                    else:
                        a = line[0] + ":" + line[1]
                    if a in chrom_pos.keys():
                        b = chrom_pos[a]
                        b = int(b) + 1
                        chrom_pos[a] = b
                    else:
                        chrom_pos[a] = 0

        #Get average number of times each variant appears
    rare_common = list(chrom_pos.values())
    print("Variant details (breed_rare_pop_common). Mean = ", mean(rare_common), "Max = ", max(rare_common), "Min = ", min(rare_common))

    with open(f"{data}/breed_rare_common_pop_chrom_pos.txt", "w") as f:
        for key in chrom_pos.keys():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)

    with open(f"Extract_variants_breed_rare_common_pop.pbs", "w") as f:
        print("#!/bin/bash -l\n"
          "#PBS -l nodes=1:ppn=8,walltime=06:00:00,mem=1g\n"
          "#PBS -m abe\n"
          "#PBS -M durwa004@umn.edu\n"
          f"#PBS -o $PBS_JOBID.Extract_variants_breed_rare_common_pop.out\n"
          f"#PBS -e $PBS_JOBID.Extract_variants_breed_rare_common_pop.err\n"
          f"#PBS -N Extract_variants_breed_rare_common_pop.pbs\n"
          "#PBS -q batch\n"
          "module load bcftools\n"
          f"cd {data}/../breed_pop_rare_common_snpeff\n"
          f"bcftools view -R {data}/breed_rare_common_pop_chrom_pos.txt ../../../joint_intersect_without_Prze/thesis_intersect.vcf.gz > breed_rare_common_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip breed_rare_common_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix breed_rare_common_pop_snpeff.vcf.gz", file = f)

    chrom_pos = {}

    for filename in os.listdir(data):
        b = filename.split(".txt")
        a = b[0]
        if "_common_pop_rare_" in b[0]:
            with open(data + "/" + filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if len(line) <2:
                        continue
                    else:
                        a = line[0] + ":" + line[1]
                    if a in chrom_pos.keys():
                        b = chrom_pos[a]
                        b = int(b) + 1
                        chrom_pos[a] = b
                    else:
                        chrom_pos[a] = 0

        #Get average number of times each variant appears
    rare_common = list(chrom_pos.values())
    print("Variant details (breed_common_pop_rare). Mean = ", mean(rare_common), "Max = ", max(rare_common), "Min = ", min(rare_common))

    with open(f"{data}/breed_common_rare_pop_chrom_pos.txt", "w") as f:
        for key in chrom_pos.keys():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)

    with open(f"Extract_variants_breed_common_rare_pop.pbs", "w") as f:
        print("#!/bin/bash -l\n"
          "#PBS -l nodes=1:ppn=8,walltime=06:00:00,mem=1g\n"
          "#PBS -m abe\n"
          "#PBS -M durwa004@umn.edu\n"
          f"#PBS -o $PBS_JOBID.Extract_variants_breed_common_rare_pop.out\n"
          f"#PBS -e $PBS_JOBID.Extract_variants_breed_common_rare_pop.err\n"
          f"#PBS -N Extract_variants_breed_common_rare_pop.pbs\n"
          "#PBS -q batch\n"
          "module load bcftools\n"
          f"cd {data}/../breed_pop_rare_common_snpeff\n"
          f"bcftools view -R {data}/breed_common_rare_pop_chrom_pos.txt ../../../joint_intersect_without_Prze/thesis_intersect.vcf.gz > breed_common_rare_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip breed_common_rare_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix breed_common_rare_pop_snpeff.vcf.gz", file = f)



    chrom_pos = {}

    for filename in os.listdir(data):
        b = filename.split(".txt")
        a = b[0]
        if "unique" in b[0]:
            with open(data + "/" + filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if len(line) <2:
                        continue
                    else:
                        a = line[0] + ":" + line[1]
                    if a in chrom_pos.keys():
                        b = chrom_pos[a]
                        b = int(b) + 1
                        chrom_pos[a] = b
                    else:
                        chrom_pos[a] = 0

        #Get average number of times each variant appears
    rare_common = list(chrom_pos.values())
    print("Variant details (unique). Mean = ", mean(rare_common), "Max = ", max(rare_common), "Min = ", min(rare_common))

    with open(f"{data}/breed_unique_chrom_pos.txt", "w") as f:
        for key in chrom_pos.keys():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)

    with open(f"Extract_variants_breed_unique_pop.pbs", "w") as f:
        print("#!/bin/bash -l\n"
          "#PBS -l nodes=1:ppn=8,walltime=06:00:00,mem=1g\n"
          "#PBS -m abe\n"
          "#PBS -M durwa004@umn.edu\n"
          f"#PBS -o $PBS_JOBID.Extract_variants_breed_unique_pop.out\n"
          f"#PBS -e $PBS_JOBID.Extract_variants_breed_unique_pop.err\n"
          f"#PBS -N Extract_variants_breed_unique_pop.pbs\n"
          "#PBS -q batch\n"
          "module load bcftools\n"
          f"cd {data}/../breed_pop_rare_common_snpeff\n"
          f"bcftools view -R {data}/breed_unique_chrom_pos.txt ../../../joint_intersect_without_Prze/thesis_intersect.vcf.gz > breed_unique_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip breed_unique_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix breed_uniqu_snpeff.vcf.gz", file = f)


