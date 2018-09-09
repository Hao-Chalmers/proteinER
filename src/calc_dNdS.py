#!/usr/bin/python

'''
This script calculates site-specific dN/dS using HyPhy's FEL one-rate output. 

HyPhy's FEL one-rate method assigns dS=1 to sites with substitutions. This script simply takes a site's dN and dS values and divides them dN/dS. This script also finds convserved sites in a sequence alignment and assigns these sites dN/dS=0.  

Author: Dariya K. Sydykova
'''

import argparse
import textwrap
import sys
from Bio import AlignIO

def calc_dNdS(aln, rates_file, out_file):

    r=open(rates_file,"r")
    out=open(out_file,"w")
    
    total_col=len(aln[0]) #total number of sites
    
    conserved_sites_lst=[]
    for i in range(total_col):
        #a list of all amino acids at site i
        col = aln[:,i]
        
        #checks if the list of amino acids at site i is identical to 
        #the first amino acid in a list repeated the number of times the list is
        if col == len(col) * col[0]:
            conserved_sites_lst.append(i+1)
        else:
            continue
    
    print(conserved_sites_lst)
    for line in r:
        if line.startswith("Site"):
            token=line.split(",")
            new_header = token[0]+",dN/dS,"+",".join(token[4:])
            out.write(new_header)
            continue
        
        token=line.split(",")
        site=token[0]
        dS=float(token[1])
        
        if dS!=0.0 or dS!=1.0:
            print("dS is not equal to 1 or 0")
            print("Check HyPhy output")
            sys.exit()
            
        dN=float(token[2])
        
        if dS==0:
            dN_dS=0
        else:
            dN_dS = dN/dS
            
        if site in conserved_sites_lst:
            new_line = site+",0,"+",".join(token[4:])
        else:
            new_line = site+","+str(dN_dS)+","+",".join(token[4:])
        
        out.write(new_line)

def main():
    '''
    Calculate site-specific dN/dS from HyPhy's FEL one-rate model output
    '''
    #creating a parser
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="Calculate site-specific dN/dS from HyPhy's FEL one-rate model output",
            epilog=textwrap.dedent('''\
            This script produces a CSV with the following columns:
            
            Column name           Description
            ===================================================================
            Site                  Site number, extracted from the alignment 
                                  FASTA file.

            dN/dS                 Site-wise dN/dS, calculated from HyPhy output.
                                  dN='beta'
                                  dS='alpha'

            LRT                   Likelihood ration test statistic for 
                                  beta = alpha, versus beta does not equal alpha
                       
            p-value               Likelihood ration test statistic for 
                                  beta = alpha, versus beta does not equal alpha 
            
            Total_branch_length   The total length of branches contributing 
                                  to inference at this site, and used to 
                                  scale dN-dS
            '''))

    #adding arguments 
    parser.add_argument('-a', metavar='<aa_aln.fasta>', type=str, help='input amino acid alignment file')
    parser.add_argument('-r', metavar='<rates.csv>', type=str, help='HyPhy FEL file')
    parser.add_argument('-o', metavar='<processed_rates.csv>', type=str, help='output processed rates file')

    args = parser.parse_args()

    #set up output file name if none is given
    if args.o is None:
        out_rates = "processed_dNdS.csv"
    else:
        out_rates = args.o
        
    aln_file=args.a
    rates_file=args.r

    aln = AlignIO.read(aln_file, "fasta") 
    calc_dNdS(aln, rates_file, out_rates)

if __name__ == "__main__":
    main()