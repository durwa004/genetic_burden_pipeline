#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:04:17 2019
@author: jonahcullen
"""

import argparse
import os
import gzip
from statistics import mean

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

    unique = {}
    rare_breed = {}
    common_breed = {}
    for filename in os.listdir(data):
        if filename.endswith(".txt"):
            b = filename.split(".txt")
            breed = b[0]
            un = filename.split("_")
            with open(data + "/" + filename, "r") as input_file:
                for line in input_file:
                    line = line.rstrip("\n").split("\t")
                    if len(line) <2:
                        continue
                    else:
                        a = line[0] + ":" + line[1]
                        if "unique" in un[0]:
                            if a in unique.keys():
                                b = unique[a]
                                b = b + ":" + breed 
                                unique[a] = b
                            else:
                                unique[a] = breed
                        rc = filename.split("common_pop")
                        if len(rc) == 2:
                            if a in common_breed.keys():
                                b = common_breed[a]
                                b = b + ":" + breed
                                common_breed[a] = b
                            else:
                                common_breed[a] = breed
                        else:
                            if a in rare_breed.keys():
                                b = rare_breed[a]
                                b = b + ":" + breed
                                rare_breed[a] = b
                            else:
                                rare_breed[a] = breed

        #Get average number of times each variant appears
    count = 0
    min_c = 3
    max_c = 0
    for key, value in common_breed.items():
        a = value.split(":")
        count += len(a)
        if len(a) < int(min_c):
            min_c = len(a)
        elif len(a) > int(max_c):
            max_c = len(a)
        
    print("Breed common pop rare variant details. Mean = ")
    print(count/(len(common_breed.keys())))
    print("Max = ")
    print(max_c)
    print("Min = ")
    print(min_c)

    count = 0
    min_c = 3
    max_c = 0
    for key, value in rare_breed.items():
        a = value.split(":")
        count += len(a)
        if len(a) < int(min_c):
            min_c = len(a)
        elif len(a) > int(max_c):
            max_c = len(a)

    print("Breed rare pop common variant details. Mean = ")
    print(count/(len(rare_breed.keys())))
    print("Max = ")
    print(max_c)
    print("Min = ")
    print(min_c)

    count = 0
    min_c = 3
    max_c = 0
    for key, value in unique.items():
        a = value.split(":")
        count += len(a)
        if len(a) < int(min_c):
            min_c = len(a)
        elif len(a) > int(max_c):
            max_c = len(a)

    print("Unique variant details. Mean = ")
    print(count/(len(unique.keys())))
    print("Max = ")
    print(max_c)
    print("Min = ")
    print(min_c)

#    with open(f"{data}/rare_common_breed_chrom_pos.txt", "w") as f, open(f"{data}/rare_common_breed_chrom_pos_with_breed_details.txt", "w") as f2:
    with open(f"{data}/all_shared_variants/breed_common_rare_pop_chrom_pos.txt", "w") as f, open(f"{data}/all_shared_variants/breed_common_rare_pop_chrom_pos_with_breed_details.txt", "w") as f2:
        for key,value in common_breed.items():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)
            print(a[0], a[1], value, sep = "\t", file = f2)

#    with open(f"{data}/rare_common_breed_chrom_pos.txt", "w") as f, open(f"{data}/rare_common_breed_chrom_pos_with_breed_details.txt", "w") as f2:
    with open(f"{data}/all_shared_variants/breed_rare_common_pop_chrom_pos.txt", "w") as f, open(f"{data}/all_shared_variants/breed_rare_common_pop_chrom_pos_with_breed_details.txt", "w") as f2:
        for key,value in rare_breed.items():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)
            print(a[0], a[1], value, sep = "\t", file = f2)


    with open(f"{data}/all_shared_variants/unique_chrom_pos.txt", "w") as f, open(f"{data}/all_shared_variants/unique_chrom_pos_with_breed_details.txt", "w") as f2:
        for key,value in unique.items():
            a = key.split(":")
            print(a[0], a[1], sep = "\t", file = f)
            print(a[0], a[1], value, sep = "\t", file = f2)

    with open(f"Extract_variants_breed_common_rare_pop.pbs", "w") as f:
        print("#!/bin/bash -l\n"
                  "#PBS -l nodes=1:ppn=8,walltime=24:00:00,mem=1g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  "#PBS -o $PBS_JOBID.Extract_variants_breed_common_rare_pop.out\n"
                  "#PBS -e $PBS_JOBID.Extract_variants_breed_common_rare_pop.err\n"
                  "#PBS -N Extract_variants_breed_common_rare_pop.pbs\n"
                  "#PBS -q batch\n"
                  "module load bcftools\n"
                  f"cd {data}/../breed_pop_rare_common_snpeff\n"
                  f"bcftools view -R {data}/all_shared_variants/breed_common_rare_pop_chrom_pos.txt ../../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz > breed_common_rare_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip breed_common_rare_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix breed_common_rare_pop_snpeff.vcf.gz", file = f)

    with open(f"Extract_variants_breed_rare_common_pop.pbs", "w") as f:
        print("#!/bin/bash -l\n"
                  "#PBS -l nodes=1:ppn=8,walltime=24:00:00,mem=1g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  "#PBS -o $PBS_JOBID.Extract_variants_breed_rare_common_pop.out\n"
                  "#PBS -e $PBS_JOBID.Extract_variants_breed_rare_common_pop.err\n"
                  "#PBS -N Extract_variants_breed_rare_common_pop.pbs\n"
                  "#PBS -q batch\n"
                  "module load bcftools\n"
                  f"cd {data}/../breed_pop_rare_common_snpeff\n"
                  f"bcftools view -R {data}/all_shared_variants/breed_rare_common_pop_chrom_pos.txt ../../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz > breed_rare_common_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip breed_rare_common_pop_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix breed_rare_common_pop_snpeff.vcf.gz", file = f)



    with open(f"Extract_variants_unique.pbs", "w") as f:
        print("#!/bin/bash -l\n"
                  "#PBS -l nodes=1:ppn=8,walltime=24:00:00,mem=1g\n"
                  "#PBS -m abe\n"
                  "#PBS -M durwa004@umn.edu\n"
                  "#PBS -o $PBS_JOBID.Extract_variants_unique.out\n"
                  "#PBS -e $PBS_JOBID.Extract_variants_unique.err\n"
                  "#PBS -N Extract_variants_unique.pbs\n"
                  "#PBS -q batch\n"
                  "module load bcftools\n"
                  f"cd {data}/../breed_pop_rare_common_snpeff\n"
                  f"bcftools view -R {data}/all_shared_variants/unique_chrom_pos.txt ../../../SnpEff/thesis_intersect_snpeff.ann.vcf.gz > unique_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip unique_snpeff.vcf && /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix unique_snpeff.vcf.gz", file = f)

