#!/usr/bin/env python3

import argparse
import os
import time
import sys

from lib import variable_alignment

def get_args():

    dirpath = os.path.dirname(__file__)
    
    parser = argparse.ArgumentParser(
        
        description='Creates a fasta of polymorphic positions, given a list of gnumbers. \
            Also outputs the number of non-variable positions per base.')
    
    parser.add_argument('-i', dest='input', required = True,
                        help='path to the input file containing one gnumber per row (no header).')
    
    parser.add_argument('-o', dest='output_prefix',required=True,
                        help='basename for output files')
    
    parser.add_argument('-s', dest='subsample', type=int, help='Subsample to N sites.')
    
    # Option to get gene-wise alignments, including non-variable sites
    #parser.add_argument('-b', dest='bed', help='Bed file with gene coordinates and names to include in the alignment')
    
    parser.add_argument('-rep', dest='repeats',
                        help='path to bed file with positions to exclude', 
                        default = os.path.join(dirpath,'resources/regions_blindspots_modlin_farhat.bed'))
    
    parser.add_argument('-dr', dest='drug',
                        help='provide the path to file containing the genomic positions you want filtered out of the fasta.',
                        default = os.path.join(dirpath,'resources/20160911_DR_filter_pos_reseqTB.txt'))
    
    parser.add_argument('-md', dest='mindepth', default=5, type = int,
                        help='Depth below which an allele is called as missing')
    
    parser.add_argument('-mm', dest='maxmissing', default=0.1, type = float,
                        help='Maximum proportion of missing alleles allowed per site')
    
    parser.add_argument('-outg',dest='outgroup', 
                        help='G number of the outgroup strain. Default is canettii ET1291 (G742339)',
                        default='G742339')
    
    parser.add_argument('-ref', dest='reference', 
                        help = 'Reference genome, used to count the number of non-variable bases.',
                        default = os.path.join(dirpath,"resources/MTB_ancestor_reference.fasta"))

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
