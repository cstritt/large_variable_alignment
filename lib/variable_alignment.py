import gzip
import numpy
import os
import struct
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
from joblib import Parallel, delayed


class mep:
    """ 
    
    """

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
        
        # Create output directory
        self.output_folder = args.output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        
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
        self.threads = args.threads


class variantmatrix:
    """ 
    
    """
    
    def __init__(self):
        
        # Filled by add_SNPs
        self.variable_positions = []
        self.variant_dict = {}
        self.outgroup_alleles = {}
        
        self.stats =  {   
            'filt_repeats' : set(),
            'filt_dr' : set(),
            'ref_larger_than_one' : 0
        }
        
        # Filled by convert_to_array
        self.mem_required = 0
        
        # Filled by traverse_depth_files
        self.n_missing = {}
        
    
    def add_SNPs(self, mep): 
        """
        Function to read in SNPs from VCF files and store them in a nested dictionary:
        
        - The first key is the position of the SNP in the MTB ancestor genome.
        - The second key is the strain name, or 'MTB_anc' for the ancestor.
        - The value is the base at that position in the given strain.
        
        """
        
        def get_variant_dict(gnumbers):
            
            variant_dict = defaultdict(dict)
            
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
                        
                        # Ignore variants in repeats and resistance loci
                        if pos in mep.repeats:
                            self.stats['filt_repeats'].add(pos)
                            continue
                        if pos in mep.dr_loci:
                            self.stats['filt_dr'].add(pos)
                            continue
                                
                        ref = fields[3]
                        alt =  fields[4].split(',')
                        
                        # If there are multiple alt alleles, take the one with the highest support
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
                            
                            variant_dict[pos][g] = alt
                                
                        # Multi-nucleotide polymorphism (optionally skip ...)
                        if len(ref) > 1 and len(ref) == len(alt):
                            self.stats['ref_larger_than_one'] += 1
                            for i, base in enumerate(alt):
                                variant_dict[pos+i][g] = base
                                
            return variant_dict
        
        self.variant_dict = get_variant_dict(mep.gnumbers)
        self.outgroup_alleles = get_variant_dict([mep.outgroup])
        
        positions = list(self.variant_dict.keys())
        self.variable_positions = sorted(positions, key=lambda x: int(x))
    
    
    def convert_to_array(self, mep):
 
        base_code = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4, 'N': 5}

        def encode_base(b):
            return base_code.get(b, 5)  # N if unknown

        samples = mep.gnumbers        
        num_sites = len(self.variable_positions)
        num_samples = len(samples)
        
        ref_bases = {pos: mep.reference[pos - 1] for pos in self.variable_positions}
        
        self.mem_required = estimate_memory(num_samples, num_sites)
        
        # Build base matrix with all "-"/0 initially, samples in rows and sites in columns
        self.matrix = numpy.full((num_samples, num_sites), base_code['-'], dtype=numpy.uint8)
        
        pos_index = {pos: i for i, pos in enumerate(self.variable_positions)}
        sample_index = {g: i for i, g in enumerate(samples)}
        
        keep_sites = []
   
        # Fill matrix
        for pos in self.variable_positions:

            i = pos_index[pos]  # column index
            ref_base = encode_base(ref_bases[pos])
            self.matrix[:,i] = ref_base
            
            samples_site = self.variant_dict[pos]
            
            for sample in samples_site:
                j = sample_index[sample]
                base = encode_base(samples_site[sample])
                self.matrix[j, i] = base
    
    
    def extract_sorted_rows_from_gz(self, path_to_depth, sorted_positions, mindepth):
  
        missing = []
        pos_idx = 0
        next_target = sorted_positions[pos_idx]
        
        with gzip.open(path_to_depth, 'rt') as f:
            for i, line in enumerate(f, start=1):
                if i == next_target:
                    val = int(line.strip())
                    if val < mindepth:
                        missing.append(i)
                    pos_idx += 1
                    if pos_idx >= len(sorted_positions):
                        break
                    next_target = sorted_positions[pos_idx]
        
        return missing
    


    def get_missing_positions(self, sample_idx, path, sorted_positions, mindepth):
        
        missing_pos = []
        
        with gzip.open(path, 'rt') as f:
            pos_idx = 0
            next_target = sorted_positions[pos_idx]
            for i, line in enumerate(f, start=1):
                if i == next_target:
                    if int(line.strip()) < mindepth:
                        missing_pos.append(next_target)
                    pos_idx += 1
                    if pos_idx >= len(sorted_positions):
                        break
                    next_target = sorted_positions[pos_idx]
                    
        return sample_idx, missing_pos
    
    def traverse_depth_files_parallel(self, mep):
        
        depth_files = [mep.filepaths[g][1] for g in mep.gnumbers if g != mep.outgroup]
        
        missing_positions = Parallel(n_jobs=mep.threads)(
            delayed(self.get_missing_positions)(
                idx, sample_path, self.variable_positions, mep.mindepth
                )
                for idx, sample_path in enumerate(depth_files)
        )

        pos_to_col = {pos:i for i, pos in enumerate(self.variable_positions)}
        
        for sample_idx, missing_pos in missing_positions:
            missing_cols = [pos_to_col[pos] for pos in missing_pos]
            self.matrix[sample_idx, missing_cols] = 4
            self.n_missing[mep.gnumbers[sample_idx]] = len(missing_cols)
            
        n, missing_outgroup = self.get_missing_positions('N', mep.filepaths[mep.outgroup][1], self.variable_positions, mep.mindepth)

        for pos in missing_outgroup:
            self.outgroup_alleles[pos][mep.outgroup] = '-'
            

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
                if row_int < mindepth:
                    missing.append(row_idx)
                    
        return missing
    
    
    def traverse_depth_files(self, mep):
        
        depth_files = {g: mep.filepaths[g][1] for g in mep.gnumbers if g != mep.outgroup}
        
        positions = self.variable_positions
        samples = mep.filepaths.keys()
        pos_index = {pos: i for i, pos in enumerate(positions)}
        sample_index = {g: i for i, g in enumerate(samples)}
        
        for g, depth_file in depth_files.items():
            j = sample_index[g]
            missing = self.extract_sorted_rows_from_gz(depth_file, positions, mep.mindepth)
            for pos in missing: 
                i = pos_index[pos]
                self.matrix[j,i] = 4
                
        
        # Same for outgroup
        missing = self.extract_rows_ra(mep.filepaths[mep.outgroup][1], positions, mep.mindepth)
        for pos in missing:
            if pos not in self.outgroup_alleles:
                self.outgroup_alleles[pos] = {}
            self.outgroup_alleles[pos][mep.outgroup] = '-'


    def max_missing_filter(self, mep):
        
        base_code = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4, 'N': 5}
        keep_sites = []
        ref_filt_count = 0
        miss_filt_count = 0
        
        num_samples = self.matrix.shape[0]
        num_sites = self.matrix.shape[1]
        
        for i in range(num_sites):
            
            site = self.matrix[:,i]
            
            # Check if only the REF is different
            valid_site = site[numpy.isin(site, [0, 1, 2, 3])]
            unique_values = numpy.unique(valid_site)
            ref_filt = True if unique_values.size == 1 else False
            
            # Remove sites with to many missing bases
            missing_prop = numpy.count_nonzero(site == base_code['-']) / num_samples
            miss_filt = True if missing_prop > mep.maxmissing else False
            
            if ref_filt or miss_filt:
                if ref_filt: 
                    ref_filt_count += 1
                if miss_filt:
                    miss_filt_count += 1
                    
            else:
                keep_sites.append(i)
                
       
        n_removed = num_sites - len(keep_sites)
        sys.stderr.write(f'{n_removed} sites did not pass filtering. Too much missing data: {miss_filt_count}. Only REF different: {ref_filt_count}\n')
        
        self.matrix = self.matrix[:, keep_sites]
        self.variable_positions = [self.variable_positions[i] for i in keep_sites]
        

    def write_output(self, mep):
        """ Output fasta with variable alignment, and txt file
        with positions that were included.
        
        # alignment
        with open(os.path.join(mep.output_folder, 'snp_alignment.fasta'), 'w') as fasta_handle:
            for g in self.seqs:
                rec = SeqRecord(Seq(self.seqs[g]), id=g, name='', description='')
                SeqIO.write(rec, fasta_handle, 'fasta')

        # positions
        with open(os.path.join(mep.output_folder, 'positions_in_alignment.tsv'), 'w') as f:
              for pos in self.sites_in_alignment:
                  f.write(str(pos) + '\n')
        
        
        """
        base_code = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4, 'N': 5}
        code_base = {v: k for k, v in base_code.items()}
        ref_bases = {pos: mep.reference[pos - 1] for pos in self.variable_positions}

        fasta_handle = open(os.path.join(mep.output_folder, 'snp_alignment.fasta'), 'w')
        pos_handle = open(os.path.join(mep.output_folder, 'positions_in_alignment.tsv'), 'w')

        for j, g in enumerate(mep.gnumbers):
            seq = ''.join(code_base[b] for b in self.matrix[j,:])
            rec = SeqRecord(Seq(seq), id=g, name='', description='')
            SeqIO.write(rec, fasta_handle, 'fasta')
       
        # Handle outgroup
        if mep.outgroup:
            outseq = []
            for pos in self.variable_positions:
                ref_base = ref_bases[pos]
                base = self.outgroup_alleles.get(pos, {}).get(mep.outgroup, ref_base)
                
                
                outseq.append(base)
            seq = ''.join(outseq)
            rec = SeqRecord(Seq(seq), id=mep.outgroup, name='', description='')
            SeqIO.write(rec, fasta_handle, 'fasta')

        # Store aligned positions
        for pos in self.variable_positions:
            pos_handle.write(str(pos) + '\n')
            
        fasta_handle.close()
        pos_handle.close()















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
        
         
    def get_seqs(self, mep, variantmatrix):
        
        # Pandas to numpy for speed
        #variantmatrix.variants = variantmatrix.variants.values
        
        self.stats = {
            'filt_maxmissing' : [],
            'biallelic' : 0,
            'multiallelic' : 0,
            'singletons' : 0,
            'invariant' : 0
        }
        
        self.stats.update(variantmatrix.stats)
        
        positions = variantmatrix.variable_positions if mep.subsample is None else variantmatrix.pos_subset
            
        for pos in positions:
            
            site = ''
            for g in mep.gnumbers:
                if g in variantmatrix.variants[pos]:
                    site += variantmatrix.variants[pos][g]
                else:
                    site += mep.reference[pos-1]
                    
            # Re-check if site is invariant
            site_bases_only = ''.join([base for base in site if base in 'ACGT'])
            if len(set(site_bases_only)) == 1:  
                self.stats['invariant'] += 1
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
            
        # Add outgroup alleles
        if mep.outgroup:
            
            self.seqs[mep.outgroup] = ''
        
            for pos in self.sites_in_alignment:
                if pos in variantmatrix.outgroup_alleles:
                    self.seqs[mep.outgroup] += variantmatrix.outgroup_alleles[pos][mep.outgroup]
                else:
                    self.seqs[mep.outgroup] += 'N'
                
                
        # Check if all seqs have same length
        seqlengths = [len(s) for s in self.seqs.values()]
        if len(set(seqlengths)) > 1:
            raise ValueError('Sequences have different lengths')
                    
                    
    def count_nonvariable(self, mep, variantmatrix):
        """ 
        Estimate the number of non-variable positions in the alignment. 
        Do not count sites excluded through subsampling or the maxmissing filter.
        """
        
        called_or_ignored = set(self.sites_in_alignment) | set(mep.repeats) | set(mep.dr_loci) | set(self.stats['filt_maxmissing']) | set(variantmatrix.subsample_excluded)
        
        for i, base in enumerate(mep.reference):
            
            pos = i+1
            
            if pos not in called_or_ignored:
                self.non_variable[base] += 1
                 
                                      
    def write_files(self, mep, get_stats=False):
 
        # alignment
        with open(os.path.join(mep.output_folder, 'snp_alignment.fasta'), 'w') as fasta_handle:
            for g in self.seqs:
                rec = SeqRecord(Seq(self.seqs[g]), id=g, name='', description='')
                SeqIO.write(rec, fasta_handle, 'fasta')

        # positions
        with open(os.path.join(mep.output_folder, 'positions_in_alignment.tsv'), 'w') as f:
              for pos in self.sites_in_alignment:
                  f.write(str(pos) + '\n')
                  
        # counts of non-variable bases                  
        with open(os.path.join(mep.output_folder, 'nonvariable_counts.tsv'), 'w') as f:
            for base in self.non_variable:
                f.write(base + '\t' + str(self.non_variable[base]) + '\n')
                
        # stats
        with open(os.path.join(mep.output_folder, 'stats.tsv'), 'w') as f:
            for k, v in self.stats.items():
                
                if not get_stats and k in ['biallelic', 'multiallelic', 'singletons']:
                    continue
                
                if k == 'filt_maxmissing':
                    row = k + '\t' + str(len(v)) + '\n'
                else:
                    row = k + '\t' + str(v) + '\n'
                    
                f.write(row)
    

def estimate_memory(nrow, ncol, dtype='int8'):
    """
    Estimate the memory requirement for a full NumPy matrix.

    Parameters:
    - nrow: Number of rows in the matrix.
    - ncol: Number of columns in the matrix.
    - dtype: Data type of the matrix elements. Default is 'int8'.

    Returns:
    - Memory requirement in bytes.
    """

    # Get the size of the specified data type in bytes
    size_of_element = numpy.dtype(dtype).itemsize

    # Calculate the total number of elements in the matrix
    num_elements = nrow * ncol

    # Calculate the total memory required
    memory_required = num_elements * size_of_element
    memory_required_gb = memory_required / (1024 ** 3)

    return round(memory_required_gb, 2)




#%% Not used

class variantmatrix_sparse:
    
    from scipy.sparse import lil_matrix

    def __init__(self, mep):
        self.stats = {   
            'filt_repeats': 0,
            'filt_dr': 0,
            'ref_larger_than_one': 0
        }

    def add_SNPs(self, mep):
        
        base_to_int = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4, "-": 5}
        
        variant_positions = set()
        g_to_idx = {g: idx for idx, g in enumerate(mep.gnumbers)}
        position_to_idx = {}
        allele_data = defaultdict(dict)

        def process_vcf(g):
            path_to_vcf = mep.filepaths[g][0]
            with open(path_to_vcf) as f:
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
                    alt = fields[4].split(',')

                    if len(alt) > 1:
                        info = fields[-1].split(':')
                        ad = [int(x) for x in info[1].split(',')]
                        ad_max = ad.index(max(ad))
                        alt = alt[ad_max - 1]
                    else:
                        alt = alt[0]

                    if len(ref) == 1 and len(alt) == 1:
                        allele_data[pos][g] = base_to_int(alt)
                        variant_positions.add(pos)
                    elif len(ref) > 1 and len(ref) == len(alt):
                        self.stats['ref_larger_than_one'] += 1
                        for i, base in enumerate(alt):
                            allele_data[pos+i][g] = base_to_int(base)
                            variant_positions.add(pos+i)

        for g in mep.gnumbers:
            process_vcf(g)

        self.variant_positions = sorted(variant_positions)
        position_to_idx = {pos: idx for idx, pos in enumerate(self.variant_positions)}
        self.position_to_idx = position_to_idx
        self.g_to_idx = g_to_idx

        # Sparse matrix initialization (LIL format for fast row-wise assignment)
        M, N = len(self.g_to_idx), len(self.position_to_idx)
        self.matrix = lil_matrix((M, N), dtype='U1')

        # Fill the sparse matrix
        for pos, alleles in allele_data.items():
            j = position_to_idx[pos]
            for g, base in alleles.items():
                i = g_to_idx[g]
                self.matrix[i, j] = base

    def get_missing_from_depth(self, mep):
        for g in mep.gnumbers:
            if g == mep.outgroup:
                continue
            path_to_depth = mep.filepaths[g][1]
            with gzip.open(path_to_depth, 'rb') as f:
                for pos in self.variant_positions:
                    offset = (pos - 1) * 5
                    f.seek(offset)
                    row = f.read(5)
                    if len(row) < 5:
                        continue
                    row_int = struct.unpack('i', row[:4])[0]
                    if row_int < mep.mindepth:
                        i = self.g_to_idx[g]
                        j = self.position_to_idx[pos]
                        self.matrix[i, j] = 5
