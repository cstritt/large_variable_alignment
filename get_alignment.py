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

    dirname = os.path.dirname(__file__) # 

    parser = argparse.ArgumentParser(
        description='Creates a fasta of polymorphic positions, given a list of gnumbers.')
    
    parser.add_argument('-i', dest='input', required = True,
                        help='path to the input file containing one gnumber per row (no header).')
    
    parser.add_argument('-o', dest='output_prefix',required=True,
                        help='basename for output files')
    
    parser.add_argument('-md', dest='mindepth', default=5, type = int,
                        help='Depth below which an allele is called as missing')
    
    parser.add_argument('-t', dest='threads', default = 1, type=int, help = 'Number of threads.')
    
    parser.add_argument('-db', dest='debug',  action = 'store_true', 
                        help = 'If set, write to stderr problematic sites to doublecheck.')
    
    parser.add_argument('-rep', dest='repeats',
                        help='path to bed file with positions to exclude', 
                        default = 'resources/regions_blindspots_modlin_farhat.bed')
    
    parser.add_argument('-dr', dest='drug',
                        help='provide the path to file containing the genomic positions you want filtered out of the fasta.',
                        default = 'resources/20160911_DR_filter_pos_reseqTB.txt')
    
    # parser.add_argument('-g',dest='gaps',
    #                     help='threshold for "-" and "N". Default is 0.9.',
    #                     default=0.9,type=float)
    
    #parser.add_argument('-outg',dest='outg',
    #                    help='provide the path to .all.var.vcf file of the outgroup to be used',
    #                    default=os.path.join(dirname,"../../../Pipeline_TB/Mcan/G00157.all.pos.vcf.gz"))
    
    #parser.add_argument('-outg-name',dest='outg_name',
    #                    help='name of the outgroup. Will be added to the fasta.',
    #                    default="Mycobacterium_canettii")
    
    #parser.add_argument('-ref', dest='reference', 
    #                    help = 'TB reference fasta',
    #                    default = "/scicore/home/gagneux/SOFT/MTB_ref_fasta/MTB_ancestor_reference.fasta")

    args = parser.parse_args()

    return args




def main():
    """
    Main function to run the large variable alignment pipeline.

    This function will load the input files (VCF, depth files, etc), 
    add the variable positions to the alignment, add missing alleles
    and write out the final alignment and stats.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    start_time = time.time()
       
    args = get_args()
    
    # Load input
    var_aln = variable_alignment(args)


    # Get variable positions
    vcf_start_time = time.time()    
    
    var_aln.add_SNPs()
    
    vcf_end_time = time.time()    
    vcf_time = vcf_end_time - vcf_start_time
    sys.stdout.write("Added variable positions in %f seconds\n" % (vcf_time))
    
    
    # Add missing alleles
    depth_start_time = time.time()
    #cpu_usage_before = psutil.cpu_percent()
    #mem_usage_before = psutil.virtual_memory().percent
    #disk_io_before = psutil.disk_io_counters()
    
    #var_aln.add_missing_bases(min_depth = 5) # parallelize this part ...
    #var_aln.process_depth_files(min_depth = 5) 
    var_aln.parallel_extract_rows(chunk_size = 200)
        
    depth_end_time = time.time()
    depth_time = depth_end_time - depth_start_time
    sys.stdout.write("Added missing alleles in %f seconds\n" % (depth_time))
    
    #cpu_usage_after = psutil.cpu_percent()
    #mem_usage_after = psutil.virtual_memory().percent
    #disk_io_after = psutil.disk_io_counters()
    #sys.stdout.write("CPU usage during missing allele addition: %f%%\n" % (cpu_usage_after - cpu_usage_before))
    #sys.stdout.write("Mem usage during missing allele addition: %f%%\n" % (mem_usage_after - mem_usage_before))
    #sys.stdout.write("i/o during missing allele addition: %f%%\n" % (disk_io_after - disk_io_before))

    # Write output
    var_aln.get_seqs()
    var_aln.write_alignment(args.output_prefix, '.')    
    var_aln.write_positions(args.output_prefix, '.')   
    var_aln.write_stats(args.output_prefix, '.')  
    
    # Write debugging info
    if args.debug:
        var_aln.write_weird_sites()
    
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    sys.stdout.write("Total time: %f seconds\n" % (elapsed_time))
    
    
if __name__ == '__main__':
    main()
