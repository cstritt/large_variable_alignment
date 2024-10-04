import gzip
import numba
import numpy as np
import os
import pandas as pd
import struct

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

class mep:
    
    def __init__(self, args):
       
        # Hard-code some parameters to v2
        self.path_to_v2_output = '/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/v2'
        self.snps_suffix = 'mutect2.filtered.homo.snps.vcf'
        self.depth_suffix = 'depth.gz'
        
        # G numbers and file paths
        self.gnumbers = []
        self.filepaths = {}
        
        with open(args.input) as f:
            for line in f:
                g =  line.strip().split('\t')[0]
                path_to_VCF =  os.path.join(self.path_to_v2_output,g[0:3],g[3:5],g[5:],g+'.{}'.format(self.snps_suffix))
                path_to_depth_file =  os.path.join(self.path_to_v2_output,g[0:3],g[3:5],g[5:],g+'.{}'.format(self.depth_suffix))
                self.gnumbers.append(g)
                self.filepaths[g] = (path_to_VCF, path_to_depth_file)
                
        # Add outgroup G number to filepaths
        self.outgroup = args.outgroup
        g = self.outgroup
        path_to_VCF =  os.path.join(self.path_to_v2_output,g[0:3],g[3:5],g[5:],g+'.{}'.format(self.snps_suffix))
        path_to_depth_file =  os.path.join(self.path_to_v2_output,g[0:3],g[3:5],g[5:],g+'.{}'.format(self.depth_suffix))
        self.filepaths[g] = (path_to_VCF, path_to_depth_file)
            
        # Reference genome
        self.reference = SeqIO.read(args.reference, 'fasta')
        self.reference = str(self.reference.seq)
        
        # Output prefix
        self.output_prefix = args.output_prefix
        
        # Positions in repeats
        self.repeats = set()
        
        with open(args.repeats,'r') as f:
            next(f)
            
            for line in f:
                fields = line.strip().split('\t')
                start = int(fields[1])
                end = int(fields[2])
                
                for pos in range(start, end + 1):
                    self.repeats.add(pos)
            
        # Positions in DR loci     
        self.dr_loci = set()

        with open(args.drug,'r') as f:
            for line in f:
                dr_pos = line.strip()
                self.dr_loci.add(int(dr_pos))
                
   
        # Filters: depth below which a site is considered missing, 
        # maximum proportion of missing alleles at a site
        self.mindepth = args.mindepth
        self.maxmissing = args.maxmissing


class variantmatrix:
    
    def __init__(self,mep):
        
        #self.variantmatrix = pd.DataFrame(columns=mep.gnumbers, index=pd.Index([], dtype=int))

        self.stats =  {   
            'filt_repeats' : 0,
            'filt_dr' : 0
        }
        

    def add_SNPs(self, mep):        
        """
        Function to read in SNPs from VCF files and store them in a nested dictionary.
        
        The nested dictionary has the following structure:
        
        - The first key is the position of the SNP in the MTB ancestor genome.
        - The second key is the strain name, or 'MTB_anc' for the ancestor.
        - The value is the base at that position in the given strain.
        
        The function also counts the number of SNPs that are multi-nucleotide 
        polymorphisms (MNPs) and stores this in the 'MNPs' key of the self.stats 
        dictionary.
        """
        
        def get_variant_dict(gnumbers):
            
            variant_dict = {}
            
            for g in gnumbers:
                
                # Get path to vcf on sciCORE
                path_to_VCF =  mep.filepaths[g][0]
        
                # Fill in variants
                with open(path_to_VCF) as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        
                        fields = line.strip().split('\t')
                        pos = int(fields[1])
                        if pos in mep.repeats:
                            self.stats['filt_repeats'] += 1
                            continue
                
                        if pos in mep.dr_loci:
                            self.stats['filt_dr'] += 1
                            continue
                                
                        ref = fields[3]
                        alt =  fields[4].split(',')
                        
                        if len(alt) > 1:
                            info = fields[-1].split(':')
                            ad = [int(x) for x in info[1].split(',')]
                            ad_max = ad.index(max(ad))
                            alt = alt[ad_max - 1]
                        else:
                            alt = alt[0]
    
                        # Single nucleotide polymorphism
                        if len(ref) == 1:
                            
                            # The ALT allele with the highest frequency is an insertion (quite rare ...): 
                            # skip, since first base of the insertion corresponds to the reference base
                            if len(alt) > 1:
                                continue
                            
                            if pos not in variant_dict:
                                variant_dict[pos] = {}     
                            variant_dict[pos][g] = alt
                                
                        # Multi-nucleotide polymorphism
                        if len(ref) > 1 and len(ref) == len(alt):

                            for i,base in enumerate(alt):
                                
                                if pos+1 not in variant_dict:
                                    variant_dict[pos+i] = {} 
                                variant_dict[pos+i][g] = base
                                
            return variant_dict
        
        self.variants = get_variant_dict(mep.gnumbers)
        self.outgroup_alleles = get_variant_dict([mep.outgroup])
        self.variable_positions = sorted(list(self.variants.keys()))
        
        # Convert nested dictionary to pandas dataframe for faster access                    
        self.variants = pd.DataFrame.from_dict(self.variants, orient='index')


    def extract_rows_ra(self, path_to_depth, positions, mindepth):
        """ Use random access to extract rows rapidly from depth file
        """
        
        with gzip.open(path_to_depth, 'rb') as f:
            
            # Get the file size
            f.seek(0, 2)
            file_size = f.tell()
            
            missing = []
            
            # Iterate over the row indices
            for row_idx in positions:
                # Calculate the offset of the row in the file
                offset = (row_idx - 1) * 5  # 4 bytes per int, +1 for newline
                
                # Check if the offset is within the file bounds
                if offset >= file_size:
                    break
                
                # Seek to the offset
                f.seek(offset)
                
                # Read the row
                row = f.read(5)
                
                # Check if we've reached the end of the file
                if len(row) < 5:
                    break
                
                # Unpack the row
                row_int = struct.unpack('i', row[:4])[0]
                
                # Add missing position to list
                if row_int > mindepth:
                    missing.append(row_idx)
                    
        return missing
    
    def traverse_depth_files(self, mep):
        
        depth_files = {g: mep.filepaths[g][1] for g in mep.gnumbers if g != mep.outgroup}
        for g, depth_file in depth_files.items():
            missing = self.extract_rows_ra(depth_file, self.variable_positions,mep.mindepth)
            self.variants.loc[missing, g] = '-' 
            
        # Same for outgroup
        missing = self.extract_rows_ra(mep.filepaths[mep.outgroup][1], self.variable_positions, mep.mindepth)
        for pos in missing:
            if pos not in self.outgroup_alleles:
                self.outgroup_alleles[pos] = {}
            self.outgroup_alleles[pos][mep.outgroup] = '-'
        self.outgroup_alleles = pd.DataFrame.from_dict(self.outgroup_alleles, orient='index')


class output:

    def __init__(self, mep):

        self.seqs = {g : '' for g in mep.gnumbers}

        self.sites_in_alignment = []
        
        self.non_variable = {
            'A' : 0,
            'C' : 0,
            'G' : 0,
            'T' : 0
        }
        
        
    def get_seqs(self, mep, variantmatrix, get_stats=False):
        
        # Pandas to numpy for speed
        variantmatrix.variants = variantmatrix.variants.values
        
        self.stats = {
            'filt_maxmissing' : [],
            'biallelic' : 0,
            'multiallelic' : 0,
            'singletons' : 0
        }
        
        self.stats.update(variantmatrix.stats)
            
        for i in range(len(variantmatrix.variants)):
            
            site_raw = variantmatrix.variants[i]
            pos = variantmatrix.variable_positions[i]
        
            site = ''

            for gt in site_raw:

                # Alt allele or missing (-)
                if isinstance(gt, str):
                    site += gt
                # Ref allele
                else:
                    site += mep.reference[pos-1]
                    
            # Skip sites for which only the reference differs
            site_no_missing = site.replace('-', '')
            if len(set(site_no_missing)) == 1:
                continue
            
            # Proportion of missing alleles
            n_missing = site.count('-')
            missing_prop = n_missing / len(mep.gnumbers)
            if missing_prop > mep.maxmissing:
                self.stats['filt_maxmissing'].append(pos)
                continue
            
            # Add bases to sequences
            for base,g in zip(site, mep.gnumbers):
                self.seqs[g] += base
                
            self.sites_in_alignment.append(pos)
            
            if get_stats:
                alleles = set(site_no_missing)
                # Allele counting                                
                if len(alleles) == 2:
                    self.stats['biallelic'] += 1
                elif len(alleles) > 2:
                    self.stats['multiallelic'] += 1
                    
                # Singleton?
                counts = Counter(site_no_missing)
                freqs = sorted(counts.values())
                if freqs == [1,len(site_no_missing)-1]:
                    self.stats['singletons'] += 1
                    
                    
        # Add outgroup alleles
        if mep.outgroup:
            
            self.seqs[mep.outgroup] = ''
        
            for pos in self.sites_in_alignment:
                if pos in variantmatrix.outgroup_alleles.index:
                    self.seqs[mep.outgroup] += variantmatrix.outgroup_alleles.loc[pos, mep.outgroup]
                else:
                    self.seqs[mep.outgroup] += 'N'
                
                
        # Check if all seqs have same length
        seqlengths = [len(s) for s in self.seqs.values()]
        if len(set(seqlengths)) > 1:
            raise ValueError('Sequences have different lengths')
                    
                    
    def count_nonvariable(self, mep):
        """ 
        Estimate the number of non-variable positions in the alignment. 
        """
        
        called_or_ignored = set(self.sites_in_alignment) | set(mep.repeats) | set(mep.dr_loci)
        
        for i, base in enumerate(mep.reference):
            
            pos = i+1
            
            if pos not in called_or_ignored:
                self.non_variable[base] += 1
                 
                                      
    def write_files(self, mep, get_stats=False):
 
        # alignment
        with open(f'{mep.output_prefix}.aligned.fasta', 'w') as fasta_handle:
                    
            for g in self.seqs:
                rec = SeqRecord(Seq(self.seqs[g]), id=g, name='', description='')
                SeqIO.write(rec, fasta_handle, 'fasta')

        # positions
        with open(f'{mep.output_prefix}.positions.tsv', 'w') as f:
              for pos in self.sites_in_alignment:
                  f.write(str(pos) + '\n')
                  
        # counts of non-variable bases                  
        with open(f'{mep.output_prefix}.nonvariable.tsv', 'w') as f:
            for base in self.non_variable:
                f.write(base + '\t' + str(self.non_variable[base]) + '\n')
                
        # stats
        with open(f'{mep.output_prefix}.stats.tsv', 'w') as f:
            for k, v in self.stats.items():
                
                if not get_stats and k in ['biallelic', 'multiallelic', 'singletons']:
                    continue
                
                if k == 'filt_maxmissing':
                    row = k + '\t' + str(len(v)) + '\n'
                else:
                    row = k + '\t' + str(v) + '\n'
                    
                f.write(row)
    