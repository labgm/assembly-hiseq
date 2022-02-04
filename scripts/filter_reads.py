#!/usr/bin/env python3

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Filters contigs, keeping only the ones above a certain threshold')
parser.add_argument('length', type=int, help='Minimum length of contigs that should be kept')
parser.add_argument('input', type=str, help='Input FASTA file')
parser.add_argument('output', type=str, help='Output FASTA file')
args = parser.parse_args()

contigs = []
for record in SeqIO.parse(args.input, "fastq"):
    if len(record.seq) > args.length:
        contigs.append(record[0:args.length]);
    elif len(record.seq)==args.length:
        contigs.append(record);

SeqIO.write(contigs, args.output, 'fastq')
