#!/usr/bin/python

'''
This script takes in protein (one for each species) and CDS sequences (multiple),
and output corresponding CDS into a new file, as well reformatted proteins in fasta

Author: Hao Wang
'''

import argparse
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import warnings

## pseudo code

# 1. load both aa and cds fasta files into SeqRecord using biopython
# 2. loop through AA sequences and 
#    get the targetSrt and search from cds description


# load both aa and cds files in fasta using biopython
faaFile = "ENSG00000117448_filtered.fasta"
fnaFile = "ENSG00000117448.fna"

# read in codon sequence file as a dictionary, key=sequence ID, value=sequence itself
# the key is same as id, while the value is a Bio.SeqRecord.SeqRecord
cds_dict = SeqIO.to_dict(SeqIO.parse(fnaFile, "fasta"))
faa_dict = SeqIO.to_dict(SeqIO.parse(faaFile, "fasta"))

# create new lists for faa and cds
new_faa_list = list()
new_cds_list = list()

# go through faa_dict
for key, value in faa_dict.items():

    # extract gene id and targetStr from description
    match = re.search(r'GeneID=(\d+)', value.description)  # the 'r' is very important
    if not match:
        value.description = ""
    else:
        newid = match.group(1)
        value.id = newid

        # extract substring for targeting cds sequence
        targetStr = value.description[value.description.find('['):]
        targetStr = targetStr.replace('isoform', 'transcript')
        value.description = ""

        # go through cds_dict and pick the corresponding one
        found = False   # see if there is any cds sequence missing
        for cdsKey, cdsValue in cds_dict.items():
            if targetStr in cdsValue.description:
                cdsValue.id = newid
                cdsValue.description = ""
                new_cds_list.append(cdsValue)
                found = True    
                break
        if not found:
            warningMsg = targetStr+" is not found in cds sequences."
            warnings.warn(warningMsg)

    new_faa_list.append(value)   # update value to new_faa
    

SeqIO.write(new_faa_list, "example.faa", "fasta")
SeqIO.write(new_cds_list, "example.fna", "fasta")
# SeqIO.write requires a list of SeqRecord objects
# if id and description are the same, use id; if not same, combine them


