#!/usr/bin/python
# coding=utf-8
"""
A script to obtain a set of unique sequences starting from sequences in FASTA format.

Usage:

    unique_seqs.py --fasta_file test_files/test_sequences.fasta --output_file test_files/unique_seqs.py

    This command will find out and write to output all the unique sequences among the dataset.
    Useful when preparing input data for phylogeny reconstruction.

"""

import argparse

# This tool needs BioPython library
from Bio import SeqIO

# Argument parser declaration
parser = argparse.ArgumentParser(prog='unique_seqs',
                                 description='Extract a subset of unique sequences, starting from a FASTA file \
                                 containing the parent sequences', )

parser.add_argument('--fasta_file', nargs=1, required=True, help='FASTA file to be processed')
parser.add_argument('--output_file', nargs=1, required=True, help='Output file')
args = parser.parse_args()

sequences_dict = SeqIO.index(args.fasta_file[0], 'fasta')
sequences = sorted(sequences_dict.keys())
print "Loaded {0:s} sequences".format(str(len(sequences)))

unique_seqs = []
duplicated_seqs = []

i = 0
j = 0

# Iterate over sequences ids
while i < len(sequences):

    # Extract master sequence for comparison
    master_seq = sequences_dict[sequences[i]]
    duplicated_seqs = []

    # Iterate over the remaining sequences, starting from the ones
    # that follow master sequence. Previous sequences are already
    # tested for uniqueness, and must not be tested again
    for j in range(i + 1, len(sequences)):

        #Extract sequence to be compared
        iter_seq = sequences_dict[sequences[j]]

        #Convert to str, otherwise comparison will fail. Also lower all
        #letters to avoid problems.
        if str(iter_seq.seq).lower() == str(master_seq.seq).lower():
            duplicated_seqs.append(iter_seq.id)
            print "Sequence {0:s} is identical to sequence {1:s}. Removing the latter.".format(master_seq.id,
                                                                                               iter_seq.id)

    # Duplicated sequences are removed *after* iteration
    # otherwise removing one may cause an IndexError while
    # iterating on j
    for seq_id in duplicated_seqs:
        sequences.remove(seq_id)

    i += 1

print "Retained {0:s} sequences".format(str(len(sequences)))
output = open(args.output_file[0], 'w+')
for seq_id in sequences:
    SeqIO.write(sequences_dict[seq_id], output, 'fasta')

output.close()