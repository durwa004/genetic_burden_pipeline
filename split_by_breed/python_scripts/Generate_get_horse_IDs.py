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
            "-i", "--ids",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="list of horse ids with breed [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    input_file = args.ids

    horse = []
    with open(input_file) as inputfile:
        inputfile.readline()
        for line in inputfile:
            line = line.rstrip("\n").split("\t")
            horse.append(line[0])

    pbs = data + "/horse_genomes_IDs.txt"
    
    with open(pbs, "w") as f:
        for i in range(len(horse)):
            print(horse[i], file = f)
