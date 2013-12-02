#!/usr/bin/python
# coding=utf-8
"""
A script to rename sequences  in FASTA format.

Usage:

    rename_fasta.py --folder test_files/Sequences/ --excel_file test_files/Tags.xlsx

    This command will rename all the sequence identifiers in the FASTA files contained in the given folder adding
    a tag, specified in the file Tags.xlsx.

    In particular, the following sequence identifier schema is supported:

    - "gi|code1|gb|code2| text" ends up in "code2"
    - "code1.x text" ends up in "code1"
    - "code1 text" ends up in "code1"

    Additionally, a sequence specific TAG is applied to the end of the new identifier, preceded by a "_"
    Tags must be specified in an Excel file, and proper columns should be specified in this script matching the ones
    containing the initial code and the assigned tag.

    If no tag is specified for a sequence, a default tag "_UNK" will be applied.

"""
import sys
import os
import argparse


#This tool needs BioPython & openpyxl libraries
from Bio import SeqIO
import openpyxl


### Snip from
### http://stackoverflow.com/a/11415816/1074046
class ReadableDir(argparse.Action):
    """
    Define a argparse action to validate a Readable Dir

    """

    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir = values[0]
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace, self.dest, prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


class ReadableFile(argparse.Action):
    """
    Define a argparse action to validate a Readable File

    """

    def __call__(self, parser, namespace, values, option_string=None):
        prospective_file = values[0]
        if not os.path.isfile(prospective_file):
            raise argparse.ArgumentTypeError("readable_file:{0} is not a valid path".format(prospective_file))
        if os.access(prospective_file, os.R_OK):
            setattr(namespace, self.dest, prospective_file)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable file".format(prospective_file))


def clean_code(dirty_code):
    """
    Given a "dirty code", extracts the GB accession number or the first code present

    :param dirty_code: string with the dirty code
    :type dirty_code: str

    :returns: A clean code
    :rtype: str
    """
    if dirty_code.startswith('gi'):
        # Matches "gi|code1|gb|code2| text"
        tokens = dirty_code.split('|')
        sequence_code = tokens[3]

        if '.' in sequence_code:
            # Remove sequence version if any
            sequence_code = sequence_code.split('.')[0]

        return sequence_code

    elif '.' in dirty_code:
        # Matches "code.x text"
        return dirty_code.split('.')[0]

    else:
        # Matches "code text"
        return dirty_code.split(' ')[0]

#
# Set some config values
#
GENBANK_CODE_COLUMN = 2
TAG_COLUMN = 5
RESULT_FOLDER_NAME = 'processed_files'
RESULT_FOLDER_NAME += '/'

# Argument parser declaration
parser = argparse.ArgumentParser(prog='rename_fasta',
                                 description='Renames FASTA identifier with only GenBank Accession Number and appends '
                                             'a tag specified in an Excel file')
parser.add_argument('--folder', nargs=1, action=ReadableDir, required=True,
                    help='Folder containing the FASTA files to be renamed')
parser.add_argument('--excel_file', nargs=1, action=ReadableFile, required=True, help='Excel file containing the tags')
args = parser.parse_args()

# Read up the Excel file
workbook = openpyxl.load_workbook(args.excel_file)
worksheet = workbook.get_sheet_names()[0]
worksheet = workbook.get_sheet_by_name(worksheet)

# Extract the columns containing the Genbank  and the assigned tag
genbank_column = worksheet.columns[GENBANK_CODE_COLUMN - 1]
tag_column = worksheet.columns[TAG_COLUMN - 1]

# Build the tags dict
tags_dict = {}
for code, tag in zip(genbank_column[1:], tag_column[1:]):
    if code.value is not None:
        tags_dict[clean_code(code.value)] = tag.value

# Go to FASTA directory and list the files
# with .fasta extension
os.chdir(args.folder)
dir_files = os.listdir('.')
fasta_files = [files for files in dir_files if files.endswith('.fasta')]

# If no FASTA files are present, exit
if len(fasta_files) == 0:
    print 'No FASTA file in given directory'
    sys.exit(0)

# Create the results dir, if missing
if not os.path.exists(RESULT_FOLDER_NAME) or not os.path.isdir(RESULT_FOLDER_NAME):
    os.makedirs(RESULT_FOLDER_NAME)

# Iterate over the files
for fasta_file in fasta_files:
    fasta = SeqIO.parse(fasta_file, 'fasta')
    new_seqs = []

    # Rename each sequence
    for sequence in fasta:
        code = clean_code(sequence.description)

        tag = ''
        try:
            tag = tags_dict[code]
        except KeyError:
            print 'Code {0:s} has no assigned tag. Defaulting to "UNK"!'.format(code)
            tag = 'UNK'

        # Clean sequence description
        sequence.description = ''

        # Assign new code, and append the sequence to the output list
        sequence.id = code + '_' + tag
        new_seqs.append(sequence)

    #Output the list of sequences in FASTA format
    output_handle = open(RESULT_FOLDER_NAME + os.path.basename(fasta_file), 'w')
    SeqIO.write(new_seqs, output_handle, 'fasta')
    output_handle.close()
