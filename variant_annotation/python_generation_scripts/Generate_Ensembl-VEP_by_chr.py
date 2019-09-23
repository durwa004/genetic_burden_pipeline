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
   
    species = "GCF_002863925.1_EquCab3.0_genomic"
    ref = "/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic.fna"

    print(f"Using {ref} and {species} for Ensembl-VEP")

    tmp = []
    for file_name in os.listdir(data):
        if file_name.endswith(".vcf.gz"):
            a = file_name.split("_intersect_joint.vcf.gz")
            tmp.append(a[0])

    for i in tmp:
        pbs = os.path.join(os.getcwd(), "Ensembl-VEP_intersect_" + i + ".pbs")
        header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=4,walltime=24:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.Ensembl-VEP_intersect_" + i + ".out\n"
              "#PBS -e $PBS_JOBID.Ensembl-VEP_intersect_" + i + ".err\n"
              "#PBS -N Ensembl-VEP_intersect_" + i + ".pbs\n"
              "#PBS -q mcqueue\n"
              "source /home/mccuem/shared/.local/conda/bin/activate ensembl-vep"
                 )

        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print("vep -i " + i + f"_intersect_joint.vcf.gz --fasta {ref} --species {species} -o Ensembl-VEP/" + i + "_Ensembl-VEP_intersect.ann --check_existing -hgvs --gene_phenotype --regulatory --plugin SpliceRegion --custom /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/variant_annotation/Ensembl-VEP/NCBI/GCF_002863925.1_EquCab3.0_genomic.modified.sorted.gff.gz,,gff,overlap,0", file=f, sep="")
