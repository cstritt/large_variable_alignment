#!/usr/bin/env python

""" To do:

    - debug mode: return histogram for missing data, per strain and per site
    -x parallelize depth file traversals
    - check indexing of different files (0 or 1)
    -x use indices rather than g numbers in dictionaries?
    - filters for max missing strain and site
    
    -x DNPs: extract correct base
    -x first evaluate site before adding bases to seq

class args:

    def __init__(self):

        self.input = '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/D_SR_data/alignment/10strains.txt'
        self.directory = '/home/cristobal/TB/projects/detettoreTB/NOTEBOOKS/D_SR_data/alignment/data'
        self.bed = '/home/cristobal/TB/gitlab/PipelineMTB/Pipeline_TB_v2/repetitive_regions/regions_blindspots_modlin_farhat.bed'
        self.drug = '/home/cristobal/TB/gitlab/PipelineMTB/Pipeline_TB/useful_scripts/20160911_DR_filter_pos_reseqTB.txt'
        self.gaps = 0.9
        self.mindepth = 5
        self.outg = ''
        self.suffix = "mutect2.filtered.homo.snps.vcf"
        self.threads = 1
        self.output_prefix = 'test'
        
args = args()

"""

import argparse
import os
import time
import sys

from lib import variable_alignment

def get_args():

    dirname = os.path.dirname(__file__)

    parser = argparse.ArgumentParser(
        
        description='Creates a fasta of polymorphic positions, given a list of gnumbers. \
            Also outputs the number of non-variable positions per base.')
    
    parser.add_argument('-i', dest='input', required = True,
                        help='path to the input file containing one gnumber per row (no header).')
    
    parser.add_argument('-o', dest='output_prefix',required=True,
                        help='basename for output files')
    
    parser.add_argument('-s', dest='subsample', type=int, help='Subsample to N strains')
    
    # Option to get gene-wise alignments, including non-variable sites
    #parser.add_argument('-b', dest='bed', help='Bed file with gene coordinates and names to include in the alignment')
    
    parser.add_argument('-rep', dest='repeats',
                        help='path to bed file with positions to exclude', 
                        default = 'resources/regions_blindspots_modlin_farhat.bed')
    
    parser.add_argument('-dr', dest='drug',
                        help='provide the path to file containing the genomic positions you want filtered out of the fasta.',
                        default = 'resources/20160911_DR_filter_pos_reseqTB.txt')
    
    parser.add_argument('-md', dest='mindepth', default=5, type = int,
                        help='Depth below which an allele is called as missing')
    
    parser.add_argument('-mm', dest='maxmissing', default=0.1, type = float,
                        help='Maximum proportion of missing alleles allowed per site')
    
    parser.add_argument('-outg',dest='outgroup', 
                        help='G number of the outgroup strain. Default is canettii ET1291 (G742339)',
                        default='G742339')
    
    parser.add_argument('-ref', dest='reference', 
                        help = 'Reference genome, used to count the number of non-variable bases.',
                        default = "resources/MTB_ancestor_reference.fasta")

    args = parser.parse_args()

    return args


def main():

    start_time = time.time()
       
    args = get_args()
    
    # Mise en place
    mep = variable_alignment.mep(args)

    # Get variable positions
    sys.stdout.write("Getting variable positions\n")
    sys.stdout.flush()
    start_varpos = time.time()
    variantmatrix = variable_alignment.variantmatrix(mep)
    variantmatrix.add_SNPs(mep)
    sys.stdout.write(f'{len(variantmatrix.variable_positions)} variable positions\n')
    sys.stdout.flush()
    end_varpos = time.time()
    sys.stdout.write("Added variable positions in %f seconds\n" % (end_varpos - start_varpos))
    sys.stdout.flush()
    
    # Add missing positions
    sys.stdout.write("Adding missing positions\n")
    sys.stdout.flush()
    start_missing = time.time()
    variantmatrix.traverse_depth_files(mep)
    end_missing = time.time()
    sys.stdout.write("Added missing positions in %f seconds\n" % (end_missing - start_missing))
    sys.stdout.flush()
    
    # Write output
    output = variable_alignment.output(mep)
    output.get_seqs(mep, variantmatrix)
    output.count_nonvariable(mep, variantmatrix)
    output.write_files(mep)
   
    # Done
    end_time = time.time()
    sys.stdout.write("Total time: %f seconds\n" % (end_time - start_time))
    sys.stdout.flush()
    
if __name__ == '__main__':
    main()
