#!/usr/bin/env python3

#%%
import argparse
import numpy as np
import sys

from Bio import SeqIO
from Bio.Seq import Seq
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
    
    parser.add_argument('-i', dest='input_fasta', required=True,
                        help='Fasta file with alignment.')
    
    parser.add_argument('-o', dest='output_fasta', required=True,
                        help='Path to filtered output fasta.')
    
    parser.add_argument('--sites', dest='site_missing', required=True,
                        help='File with the count/proportion of missing alleles per site.')
    
    parser.add_argument('--samples', dest='sample_missing', required=True,
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
    outhandle = open(args.output_fasta, 'w')
    
    # Remove samples with too many missing sites
    samples_to_remove = set()
    with open(args.sample_missing) as f:
        for line in f:
            fields = line.strip().split('\t')
            sample = fields[0]
            prop_missing = float(fields[-1])
            if prop_missing > args.maxmissing_sample:
                samples_to_remove.add(sample)
    sys.stderr.write(f'Removing {len(samples_to_remove)} samples below threshold {args.maxmissing_sample}\n')
    
    # Remove sites with too many missing alleles
    sites_to_remove = set()
    with open(args.site_missing) as f:
        indx = 0
        for line in f:
            fields = line.strip().split('\t')
            sample = fields[0]
            prop_missing = float(fields[-1])
            if prop_missing > args.maxmissing_site:
                sites_to_remove.add(indx)
            indx += 1
    sys.stderr.write(f'Removing {len(sites_to_remove)} sites below threshold {args.maxmissing_site}\n')
    
    # Remove sites
    seq_length = len(next(iter(aln_d.values())).seq)
    keep_indices = np.array([i for i in range(seq_length) if i not in sites_to_remove])
    for sample, record in aln_d.items():
        if sample in samples_to_remove:
            continue
        
        seq_arr = np.array(list(record.seq))
        filtered_seq = ''.join(seq_arr[keep_indices])
        record.seq = Seq(filtered_seq)
        SeqIO.write(record, outhandle, 'fasta')
        
    outhandle.close()
        
if __name__ == '__main__':
    main()
