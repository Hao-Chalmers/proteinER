#!/usr/bin/python

'''
This script extracts dN/dS values from HyPhy output file for the first sequence in the alignment

Author: Dariya K. Sydykova
'''

import argparse
from Bio import AlignIO

def extract_dNdS(aln_file, rates_file, out_rates):	
	aln = AlignIO.read(aln_file, "fasta") #read fasta file
	first_seq=aln[0].seq #get the first sequence in the alignment
	
	r=open(rates_file,"r")
	out=open(out_rates,"w")
	
	site=0 #set up a counter to parse the alignment
	pos=1 #set up a counter for sequence position
	for line in r:
		if line.startswith("dN/dS"):
			out.write('fasta_position\tfasta_aa'+line) #write a new heading
			continue
		
		if first_seq[site]!='-': #if the site is not a gap, write the fasta position, amino acid, and dN/dS value to the output file
			aa=first_seq[site]
			new_line=str(pos)+','+aa+','+line
			out.write(new_line)
			pos+=1
			
		site+=1		

def main():
	'''
	Extract dN/dS values from HyPhy output for non-gap sites
	'''
	#creating a parser
	parser = argparse.ArgumentParser(description='Assign dN/dS of 0 to conserved sites.')
	#adding arguments 
	parser.add_argument('-a', metavar='<aa_aln.fasta>', type=str, help='input amino acid alignment file')
	parser.add_argument('-r', metavar='<rates.csv>', type=str, help='HyPhy FEL1 file')
	parser.add_argument('-o', metavar='<processed_rates.csv>', type=str, help='output processed rates file')

	args = parser.parse_args()

	#set up output file name if none is given
	#set up output file name if none is given
	if args.o is None:
		out_rates = "processed_"+args.r
	else:
		out_rates = args.o
		
	aln_file=args.a
	rates_file=args.r

	extract_dNdS(aln_file, rates_file, out_rates)
	

if __name__ == "__main__":
	main()