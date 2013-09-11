#!/usr/bin/python
# coding=utf-8
"""
A script to extract identical subsequences starting from sequences in FASTA format.

Usage:

    extract_subsequence.py --fasta_file test_files/test_sequences.fasta 10 20

    This command will extract the range nucleotides from position 10 to 20 from all the sequences
    contained in the fasta file test_sequences.fasta .


    If you need to extract a subsequence from a given start position to the end of the seq, simply specify only
    the starting point:

    extract_subsequence.py --fasta_file test_files/test_sequences.fasta 10


    You can combine more than a subsequence in a single, composite, subsequence:

    extract_subsequence.py --fasta_file test_files/test_sequences.fasta 10 20 30 40

    this will output a single subsequence composed by the subsequences 10-20 30-40 for each sequence in the FASTA file.


    If end position is smaller than start position, the script will assume that you are trying to extract a subsequence
    that encompasses the end of the mtDNA molecule (e.g. Dloop). In that case, the program will combine two "virtual"
    subsequences:

            start_point-endofsequence & startofsequence-end_point

    Obviously, also in this case more than a subsequence can be extracted at a time.


"""

import argparse

#This tool needs BioPython library
from Bio import SeqIO, SeqRecord, Seq


def extractAroundEnd(sequence, startPosition, endPosition):
    """
    Extract a subsequence that encompasses the sequence end.


    :param sequence: Parent sequence
    :type sequence: SeqObject
    :param startPosition: Subsequence start position (should be greater than endPosition)
    :type startPosition: int
    :param endPosition: Subsequence end position (should be smaller than startPosition
    :type endPosition: int

    :returns: Subsequence of parent sequence
    :rtype: SeqObject

    """

    return sequence[startPosition - 1:] + sequence[:endPosition]


def extractSubsequence(sequence, startPosition, endPosition=None):
    """
    Extract a subsequence.


    :param sequence: Parent sequence
    :type sequence: SeqObject
    :param startPosition: Subsequence start position (must be smaller than endPosition)
    :type startPosition: int
    :param endPosition: Subsequence end position (must be smaller than startPosition
    :type endPosition: int

    :returns: Subsequence of parent sequence
    :rtype: SeqObject

    """
    if endPosition:
        return sequence[startPosition - 1:endPosition]
    else:
        return sequence[startPosition - 1:]


# Argument parser declaration
parser = argparse.ArgumentParser(prog='extract_subsequence',
                                 description='Extract a subsequence from a set of sequences, eventually combining more \
                                 than one subsequence',
)

parser.add_argument('--fasta_file', nargs=1, type=file, required=True, help='FASTA file to be processed')
parser.add_argument('--output_file', nargs=1, required=True, help='Output file')
parser.add_argument('subsequences', nargs='+', type=int, help='Subsequences to be extracted and joined. Eg: 10 20')
args = parser.parse_args()

# Open handle for output
output = open(args.output_file[0], 'w+')

sequences = SeqIO.parse(args.fasta_file[0], 'fasta')
subsequences_length = len(args.subsequences)

# If len requestes subsequences is odd, the last subsequence must
# go to the end of the seq
if len(args.subsequences) % 2 != 0:
    last_subseq_start = args.subsequences.pop(len(args.subsequences) - 1)
else:
    last_subseq_start = None

for seq in sequences:
    subseqs = []

    for i in xrange(0, len(args.subsequences), 2):
        start_position, end_position = args.subsequences[i], args.subsequences[i + 1]

        if end_position < start_position:
            subseqs.append(extractAroundEnd(seq.seq, start_position, end_position))
        else:
            subseqs.append(extractSubsequence(seq.seq, start_position, end_position))

    if last_subseq_start:
        subseqs.append(extractSubsequence(seq.seq, last_subseq_start))

    merged_seq = SeqRecord.SeqRecord(Seq.Seq(''), id=seq.id, description=seq.description)

    for subseq in subseqs:
        merged_seq += subseq

    SeqIO.write(merged_seq, output, 'fasta')

args.fasta_file[0].close()
output.close()
