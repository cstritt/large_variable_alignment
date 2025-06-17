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

    def __init__(self, args):
        
        # Create list of vcf_file(s)
        self.vcf_files = []
        with open(args.input) as f:
            for line in f:
                self.vcf_files.append(line.strip())
                
        # Depth files, if provided
        self.depth_files = {}
        if args.depth_files:
            with open(args.depth_files) as f:
                for line in f:
                    sample, filepath = line.strip().split('\t')
                    self.depth_files[sample] = filepath
            
        # Outgroup vcf and depth file
        self.outvcf = args.outvcf
        self.outdepth = args.outdepth
            
        # Reference genome
        self.reference = SeqIO.read(args.reference, 'fasta')
        self.reference = str(self.reference.seq)
        
        # Create output directory
        self.output_folder = args.output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        
        # Positions in repeats
        self.exclude = set()
        if args.exclude:
            with open(args.exclude,'r') as f:
                next(f)
                for line in f:
                    fields = line.strip().split('\t')
                    start = int(fields[1])
                    end = int(fields[2])
                    for pos in range(start, end + 1):
                        self.exclude.add(pos)        
   
        # Filters: depth below which a site is considered missing, 
        # maximum proportion of missing alleles at a site
        self.mindepth = args.mindepth
        self.maxmissing = args.maxmissing
        self.threads = args.threads

    
class variantmatrix:
    
    def __init__(self):
        
        # Filled by add_SNPs
        self.variable_positions = []
        self.variant_dict = {}
        self.outgroup_alleles = {}
        
        self.stats =  {  
            'SNPs': 0,
            'MNPs': 0,
            'insertions' : 0,
            'deletions' : 0,
            'exclude' : set()
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
        
        def get_variant_dict(vcf_files):
            
            # Create a dictionary to store the variants
            variant_dict = defaultdict(dict)
            samples_all = []
                       
            # Loop through the vcf files
            for vcf in vcf_files:
                
                if vcf.endswith('.gz'):
                    vcfhandle = gzip.open(vcf, 'rt')
                else:
                    vcfhandle = open(vcf)
                    
                for line in vcfhandle:
                    
                    if line.startswith('##'):
                        continue
                    
                    fields = line.strip().split('\t')
                    
                    if line.startswith('#CHROM'):
                        samples = fields[9:]
                        samples_all += samples
                        continue
                    
                    pos = int(fields[1])
                    
                    # Ignore variants in excluded loci
                    if pos in mep.exclude:
                        self.stats['exclude'].add(pos)
                        continue
    
                    # Get the REF and ALT alleles
                    ref = fields[3]
                    alt =  fields[4].split(',')
                    
                    # If there are multiple alt alleles, take the one with the highest support
                    # Presupposes that vcf are single-sample or have been merged with each allele on a seperate row
                    if len(alt) > 1:
                        info = fields[-1].split(':')
                        ad = [int(x) for x in info[1].split(',')]
                        ad_max = ad.index(max(ad))
                        alt = alt[ad_max - 1]
                    else:
                        alt = alt[0]
                                
                    genotypes = [x.split(':')[0] for x in fields[9:]]
                    genotypes_valid = [i for i, gt in enumerate(genotypes) if not '.' in gt]
                    
                    for i in genotypes_valid:
                        
                        gt = genotypes[i]
                        sample = samples[i]
                        
                        # Insertion
                        if len(ref) == 1 and len(alt) > 1:
                            self.stats['insertions'] += 1
                            continue
                        
                        # SNP
                        if len(ref) == 1 and len(alt) == 1:
                            variant_dict[pos][sample] = alt
                            self.stats['SNPs'] += 1
                            
                        # MNP
                        if len(ref) > 1 and len(ref) == len(alt):
                            for i, base in enumerate(alt):
                                # Add the variant to the dictionary
                                variant_dict[pos+i][sample] = base
                            self.stats['MNPs'] += 1
                        
                        # Deletion (redundant when depth files are provided: '-' will be overwritten with 'N')
                        if len(ref) > 1 and len(alt) == 1:
                            for i in range(1,len(ref)):
                                variant_dict[pos+i][sample] = '-'
                            
                            self.stats['deletions'] += 1
                                    
            return variant_dict, samples_all
        
        self.variant_dict, self.samples = get_variant_dict(mep.vcf_files)  # this defines the samples order for the variant matrix!
        self.outgroup_alleles, self.outgroup = get_variant_dict([mep.outvcf])
        self.outgroup = self.outgroup[0]
        
        positions = list(self.variant_dict.keys())
        self.variable_positions = sorted(positions, key=lambda x: int(x))
    
    
    def convert_to_array(self, mep):
 
        base_code = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4, 'N': 5}

        def encode_base(b):
            return base_code.get(b, 5)  # N if unknown
     
        num_sites = len(self.variable_positions)
        num_samples = len(self.samples)
        
        ref_bases = {pos: mep.reference[pos - 1] for pos in self.variable_positions}
        
        self.mem_required = estimate_memory(num_samples, num_sites)
        
        # Build base matrix with all "-"/0 initially, samples in rows and sites in columns
        self.matrix = numpy.full((num_samples, num_sites), base_code['-'], dtype=numpy.uint8)
        
        pos_index = {pos: i for i, pos in enumerate(self.variable_positions)}
        sample_index = {g: i for i, g in enumerate(self.samples)}
           
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
        
        depth_files = [mep.depth_files[sample] for sample in self.samples]
         
        missing_positions = Parallel(n_jobs=mep.threads)(
            delayed(self.get_missing_positions)(
                idx, sample_path, self.variable_positions, mep.mindepth
                )
                for idx, sample_path in enumerate(depth_files)
        )

        pos_to_col = {pos:i for i, pos in enumerate(self.variable_positions)}
        
        for sample_idx, missing_pos in missing_positions:
            missing_cols = [pos_to_col[pos] for pos in missing_pos]
            self.matrix[sample_idx, missing_cols] = 5
            self.n_missing[self.samples[sample_idx]] = len(missing_cols)
            
        n, missing_outgroup = self.get_missing_positions('N', mep.outdepth, self.variable_positions, mep.mindepth)

        for pos in missing_outgroup:
            self.outgroup_alleles[pos][self.outgroup] = 'N'
            
    
    def count_missing(self):
        """
        Count the occurrences of 4s (-) and 5s (N) in each row and column of the variant matrix.
        
        self.matrix: variant matrix with samples as rows and sites as columns.

        Returns:
        - dict: A dictionary with 'row_counts' and 'col_counts' as keys.
                'row_counts' is a list of tuples (count_4, count_5) for each row.
                'col_counts' is a list of tuples (count_4, count_5) for each column.
        """
        # Count 4s and 5s in each row
        self.sample_missing = [(numpy.sum(row == 4), numpy.sum(row == 5)) for row in self.matrix]

        # Count 4s and 5s in each column
        self.site_missing = [(numpy.sum(col == 4), numpy.sum(col == 5)) for col in self.matrix.T]         
            

    def apply_filters(self, mep):
        
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
        """ Output fasta with variable alignment, a sites file with the number of positions included 
        and the number of deleted/missing alleles per site, and a sample file with the number of missing 
        alleles per sample.
        
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
        sample_handle = open(os.path.join(mep.output_folder, 'snp_alignment.samples.tsv'), 'w')
        site_handle = open(os.path.join(mep.output_folder, 'snp_alignment.sites.tsv'), 'w')

        for j, g in enumerate(self.samples):
            seq = ''.join(code_base[b] for b in self.matrix[j,:])
            rec = SeqRecord(Seq(seq), id=g, name='', description='')
            SeqIO.write(rec, fasta_handle, 'fasta')
       
        # Handle outgroup
        if self.outgroup:
            outseq = []
            for pos in self.variable_positions:
                ref_base = ref_bases[pos]
                base = self.outgroup_alleles.get(pos, {}).get(self.outgroup, ref_base)
                outseq.append(base)
            seq = ''.join(outseq)
            rec = SeqRecord(Seq(seq), id=self.outgroup, name='', description='')
            SeqIO.write(rec, fasta_handle, 'fasta')
        fasta_handle.close()

        # Store aligned positions
        n_sites = len(self.variable_positions)
        n_samples = len(self.samples)
        
        for pos, count_missing in zip(self.variable_positions, self.site_missing):
            site_handle.write(f'{pos}\t{count_missing}\t{round(count_missing/n_sites)}\n')
        site_handle.close()
            
        for sample, count_missing in zip(self.samples, self.sample_missing):
            sample_handle.write(f'{sample}\t{count_missing}\t{round(count_missing/n_samples)}')
        sample_handle.close()
    

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
