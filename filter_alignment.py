#!/usr/bin/env python3

#%%
import argparse
import sys

from Bio import SeqIO
from collections import defaultdict
from itertools import zip_longest


class GETARGS:
    def __init__(self):
        self.input_fasta = '/scicore/home/gagneux/stritt0001/TB/projects/MTBC_IS6110-poly/workflow/results/10k/alignment/snp_alignment.fasta'
        self.output_dir = '/scicore/home/gagneux/stritt0001/github/large_variable_alignment/testing/'
        self.maxmissing_indv = 0
        self.maxmissing_site = 0
    
def get_args():
    
    parser = argparse.ArgumentParser(
        description = """
        Get distributions of missing alleles per site and sample and/or apply filters to alignment. 
        """)
    
    parser.add_argument('-i', dest='input_fastsa', required=True,
                        help='Fasta file with alignment.')
    
    parser.add_argument('-i', dest='site_missing', required=True,
                        help='File with the count/proportion of missing alleles per site.')
    
    parser.add_argument('-i', dest='sample_missing', required=True,
                        help='File with the count/proportion of missing alleles per sample.')
        
    parser.add_argument('-mm', dest='maxmissing_site', required=True, type = float,
                        help='Maximum proportion of missing alleles per site')
    
    parser.add_argument('-ms', dest='maxmissing_sample', required=True, type = float,
                        help='Maximum proportion of missing alleles per sample')
    
    args = parser.parse_args()

    return args


def main():
    
    args = get_args()
    aln_d = SeqIO.to_dict(SeqIO.parse(args.input_fasta, "fasta"))

    sites_to_remove = []
    with open(args.site_missing) as f:
        indx = 0
        for line in f:
            fields = line.strip().split('\t')
            sample = fields[0]
            prop_missing = float(fields[-1])
            if prop_missing > args.maxmissing_site:
                sites_to_remove.append(indx)
            indx += 1
    
    samples_to_remove = []
    with open(args.sample_missing) as f:
        for line in f:
            fields = line.strip().split('\t')
            sample = fields[0]
            prop_missing = float(fields[-1])
            if prop_missing > args.maxmissing_sample:
                samples_to_remove.append(sample)
                
    # Remove samples
    aln_d_filt = {k : aln_d[k] for k in aln_d if k not in samples_to_remove}
    
    # Remove sites
    for sample, record in aln_d_filt.items():
        seq = record.seq
        filtered_seq = ''.join([seq[i] for i in range(len(seq)) if i not in sites_to_remove])
        aln_d_filt[sample].seq = filtered_seq
        
    # Write to stdout
    for sample, record in aln_d_filt.items():
        SeqIO.write(record, sys.stdout, 'fasta')

        
if __name__ == '__main__':
    main()
