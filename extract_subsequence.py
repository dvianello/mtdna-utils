#!/usr/bin/python
__author__ = 'Dario Vianello'

import sys
import argparse

#This tool need BioPython library
from Bio import SeqIO


def extractDloop(sequence, startPosition, endPosition):
    return sequence[startPosition - 1:] + sequence[:endPosition]


def extractSubsequence(sequence, startPosition, endPosition):
    return sequence[startPosition - 1:endPosition]


parser = argparse.ArgumentParser(prog="extract_subsequence",
                                 description="Extract a subsequence from a set of sequences, eventually combining more \
                                 than one subsequence")

parser.add_argument('--fasta_file', nargs=1, type=file, required=True, help="FASTA file to be processed")
parser.add_argument('subsequences', nargs='+', help="Subsequences to be extracted and joined")
args = parser.parse_args()
print args

raise SystemExit

f = open(sys.argv[1], 'r')
output = open('dloops.fasta', 'w')
sequences = SeqIO.parse(f, 'fasta')
for items in sequences:
    print >> output, '>' + items.id
    if endPosition < startPosition:
        print >> output, extractDloop(items.seq, startPosition, endPosition)
    else:
        print >> output, extractSubsequence(items.seq, startPosition, endPosition)

f.close()
output.close()