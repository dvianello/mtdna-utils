#!/usr/bin/python
__author__ = 'Dario Vianello'

import sys

from Bio import SeqIO


def extractDloop(sequence, startPosition, endPosition):
    return sequence[startPosition - 1:] + sequence[:endPosition]


def extractSubsequence(sequence, startPosition, endPosition):
    return sequence[startPosition - 1:endPosition]


if len(sys.argv) != 4:
    print "Usage: extract_subsequence.py <fasta> <start> <end>"
    sys.exit(1)

inputFile = sys.argv[1]
startPosition = int(sys.argv[2])
endPosition = int(sys.argv[3])

#endPosition = int(sys.argv[2])
#startPosition = int(sys.argv[3])-1

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