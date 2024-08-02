#!/usr/bin/env python

import bioinfo
import matplotlib.pyplot as plt
import gzip
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Setting up fastq and read length arguments")
    parser.add_argument("-f", help="input fastq file", type=str)
    parser.add_argument("-l", help="sequence length", type=int)
    parser.add_argument("-r", help="read label", type=str)
    return parser.parse_args()
args = get_args()

read = []

for i in range(args.l):
    read.append(0.0)

with gzip.open(args.f, "rt") as fq:
    for line_num, line in enumerate(fq):
        line = line.strip()
        if line_num % 4 == 3:
            for base_pos, base in enumerate(line):
                read[base_pos] += bioinfo.convert_phred(base)
                
#calculate and print mean value at each base position
for index, value in enumerate(read):
    read[index] = value / ((line_num + 1) / 4)
    print(f"{index} {read[index]}")

filename = args.f.split("/")[-1].split(".")[0]

plt.title("Read"+args.r+" Per Base Average Quality Score Distribution")
plt.scatter(range(args.l), read)
plt.xlabel("Base Position")
plt.ylabel("Quality Score")
plt.savefig(filename+"_distributionplot.png")