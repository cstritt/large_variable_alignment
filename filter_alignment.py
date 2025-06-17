#!/usr/bin/env python3

#%%
import argparse

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
                        help='.')
    
    parser.add_argument('-o', dest='output_fasta',required=True,
                        help='.')
    
    parser.add_argument('-m', dest='missing_stats', )
        
    parser.add_argument('-mi', dest='maxmissing_indv', default=1, type = float,
                        help='')
    
    parser.add_argument('-ms', dest='maxmissing_site', default=0.5, type = float,
                        help='Maximum proportion of missing alleles allowed per site')
    
    args = parser.parse_args()

    return args


def count_N_and_minus_per_sample_and_position(alignment):
    """
    Count the occurrences of 'N' and '-' in nucleotide sequences from an alignment,
    both per sample and per position in the alignment.

    Args:
    - alignment (dict): A dictionary containing the alignment. Keys are sample headers and values are sequences.

    Returns:
    - dict: A dictionary with sequence headers as keys and a tuple (count_N, count_minus) as values.
    - dict: A dictionary with position indices as keys and a tuple (count_N, count_minus) as values.
    """
    sample_counts = {}
    position_counts = defaultdict(lambda: [0, 0])  # [count_N, count_minus]
    
    # Count per sample
    for header, record in alignment.items():
        sequence = record.seq
        count_N = sequence.count('N')
        count_minus = sequence.count('-')
        sample_counts[header] = (count_N, count_minus)
        
    # Count per position using zip_longest to handle sequences of different lengths
    for position in zip_longest(*alignment.values(), fillvalue=' '):
        count_N = sum(1 for nucleotide in position if nucleotide == 'N')
        count_minus = sum(1 for nucleotide in position if nucleotide == '-')
        idx = position_counts[next(iter(position_counts))][0]  # Use the next available index
        position_counts[idx] = [count_N, count_minus]
    
    return sample_counts, dict(position_counts)


def site_filter(alignment, threshold=0.1):
    """
    Filter out sites and samples where the proportion of '-' or 'N' characters exceeds a threshold.

    Args:
    - alignment (dict): A dictionary containing the alignment. Keys are sample headers and values are sequences.
    - threshold (float): The threshold proportion (between 0 and 1).

    Returns:
    - dict: Filtered alignment.
    """
    filtered_alignment = {}
    positions_to_remove = set()

    # Determine positions to remove
    for seq in alignment.values():
        for idx, nucleotide in enumerate(seq):
            if nucleotide in ['-', 'N']:
                if seq.count(nucleotide, idx, idx+1) / len(seq) > threshold:
                    positions_to_remove.add(idx)

    # Filter sequences
    for header, seq in alignment.items():
        filtered_seq = ''.join([seq[i] for i in range(len(seq)) if i not in positions_to_remove])
        if filtered_seq:  # Ensure the sample isn't empty after filtering
            filtered_alignment[header] = filtered_seq

    return filtered_alignment


def remove_samples_with_threshold(alignment, threshold=0.1):
    """
    Remove samples where the proportion of '-' or 'N' characters exceeds a threshold.

    Args:
    - alignment (dict): A dictionary containing the alignment. Keys are sample headers and values are sequences.
    - threshold (float): The threshold proportion (between 0 and 1).

    Returns:
    - dict: Filtered alignment.
    """
    filtered_alignment = {}

    for header, seq in alignment.items():
        count_N = seq.count('N')
        count_minus = seq.count('-')
        total = len(seq)
        
        if (count_N + count_minus) / total <= threshold:
            filtered_alignment[header] = seq

    return filtered_alignment

# Example usage:
# sample_results, position_results = count_N_and_minus_per_sample_and_position("path_to_your_fasta_file.fasta")
# filtered_alignment = filter_alignment_by_threshold(sample_results, threshold=0.2)
# filtered_samples = remove_samples_with_threshold(sample_results, threshold=0.2)

#%%
def main():
    
    #%%
    #args = get_args()
    args = GETARGS()
    aln_d = SeqIO.to_dict(SeqIO.parse(args.input_fasta, "fasta"))
    
    #%% Count missing and deleted alleles per site and sample
    count_N_and_minus_per_sample_and_position(aln_d)
    
    
    
    
    #%% Filter: first sites, then samples
    aln_site_filt = site_filter()
    aln_sample_filt = remove_samples_with_threshold()
        
        
if __name__ == '__main__':
    main()
