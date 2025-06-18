#!/usr/bin/env python3

import argparse
import os
import time
import sys

from lib import variable_alignment

def get_args():

    dirpath = os.path.dirname(__file__)
    
    parser = argparse.ArgumentParser(
        description = """
        Creates a fasta of polymorphic positions and a list of these positions, 
        given a list of gnumbers. Also outputs the number of non-variable positions 
        per base.
        """)
    
    parser.add_argument('-i', dest='input', required=True,
                        help='A list of one multi-sample vcf or multiple single-sample vcfs.')
    
    parser.add_argument('-o', dest='output_folder',required=True,
                        help='Path to output folder.')
        
    parser.add_argument('-e', dest='exclude',
                        help='Path to bed file with positions to exclude', 
                        default = '')
    
    parser.add_argument('-d', dest='depth_files',
                        help='Tab separated file with sample name and path to file indicating sequencing depth at each position.')
    
    parser.add_argument('-ov',dest='outvcf', 
                        help='Path to vcf containing outgroup variants. Default is canettii ET1291 (G742339)',
                        default=os.path.join(dirpath, 'resources/outgroup/G742339.mutect2.filtered.allvariants.AF_09.vcf.gz'))
    
    parser.add_argument('-od',dest='outdepth', 
                        help='Path to depth file for outgroup. Default is canettii ET1291 (G742339)',
                        default=os.path.join(dirpath, 'resources/outgroup/G742339.depth.gz'))
    
    parser.add_argument('-ref', dest='reference', 
                        help = 'Reference genome, used to count the number of non-variable bases.',
                        default = os.path.join(dirpath,"resources/reference/MTB_ancestor_reference.fasta"))
    
    parser.add_argument('-t', dest='threads', default=8, type=int,
                        help='Number of threads for parallel traversal of depth files')
    
    parser.add_argument('-md', dest='mindepth', default=5, type = int,
                        help='Depth below which an allele is called as missing')
    
    args = parser.parse_args()

    return args


def main():

    start_time = time.time()
       
    args = get_args()
    
    # Mise en place
    mep = variable_alignment.mep(args)

    # Get variable positions
    sys.stderr.write("Getting variable positions\n")
    sys.stderr.flush()
    start_varpos = time.time()
    variantmatrix = variable_alignment.variantmatrix()
    variantmatrix.add_SNPs(mep)
    sys.stderr.write(f'{len(variantmatrix.variable_positions)} variable positions\n')
    sys.stderr.flush()
    end_varpos = time.time()
    sys.stderr.write("Added variable positions in %f seconds\n" % (end_varpos - start_varpos))
    sys.stderr.flush()
    
    # Convert to numpy array and remove REF only sites
    variantmatrix.convert_to_array(mep)
    sys.stderr.write(f'Estimated memory requirement for {len(variantmatrix.samples)}x{len(variantmatrix.variable_positions)} matrix: {variantmatrix.mem_required}G\n')
    
    # Add missing positions
    if mep.depth_files:
        sys.stderr.write("Adding missing positions\n")
        sys.stderr.flush()
        start_missing = time.time()
        variantmatrix.traverse_depth_files_parallel(mep)
        end_missing = time.time()
        sys.stderr.write("Added missing positions in %f seconds\n" % (end_missing - start_missing))
        sys.stderr.flush()
    
    # Count missing alleles per site and sample
    variantmatrix.count_missing()
 
    # Write output
    variantmatrix.write_output(mep)
    
    # Done
    end_time = time.time()
    sys.stderr.write("Total time: %f seconds\n" % (end_time - start_time))
    sys.stderr.flush()
    
if __name__ == '__main__':
    main()
