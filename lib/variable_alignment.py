import concurrent.futures
import gzip
import os
import struct
import sys

class variable_positions:
    
    def __init__(self):
        pass
        
    def default_params(self):

        self.reference = 'MTB_anc'
        self.reference_length = 4411532
        self.path_to_v2_output = '/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/v2'
        self.snps_suffix = 'mutect2.filtered.homo.snps.vcf'
    
    

        
        ## Other parameters
        self.nthreads = args.threads
        
        ## G numberspass
        self.gnumbers = []
        with open(args.input) as f:
            for line in f:
                g =  line.strip().split('\t')[0]
                self.gnumbers.append(g)
        
        
        ## Positions in repeats
        self.repeats = set()
        
        with open(args.repeats,'r') as f:
            next(f)
            
            for line in f:
                fields = line.strip().split('\t')
                start = int(fields[1])
                end = int(fields[2])
                
                for pos in range(start, end + 1):
                    self.repeats.add(pos)
            

        ## Positions in DR loci     
        self.dr_loci = set()

        with open(args.drug,'r') as f:
            for line in f:
                dr_pos = line.strip()
                self.dr_loci.add(int(dr_pos))
                
   
        ## Depth below which a base is considered missing
        self.minDepth = int(args.mindepth) 
        
        
        # OUTPUT                
        self.sites_in_alignment = []
        
        
        self.stats = {
            
            #'singletons' : 0,
            'biallelic' : 0,
            'multiallelic' : 0,
            
            # filtered positions
            'filt_repetitive' : 0,
            'filt_dr_loci' : 0,
            'filt_site_missing' : 0,
            'filt_strain_missing': 0,
            
            # multi-nucleotide polymorphisms
            'MNPs' : 0
            }
        
        self.sites_to_doublecheck = {
            'Number of alleles different from number of strains' :[],
            'ALT allele and missing' : [],
            'missing' : []
        }
        
        # self.strain_stats = {}
        

    def add_SNPs(self):        
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
        
        self.variants = {}

        
        for g in self.gnumbers:
        
            # Get path to vcf on sciCORE
            path_to_VCF =  os.path.join(self.path_to_v2_output,g[0:3],g[3:5],g[5:],g+'.{}'.format(self.snps_suffix))
            #path_to_VCF =  os.path.join(self.path_to_v2_output,g+'.{}'.format(self.snps_suffix))
            
            # Fill in variants
            with open(path_to_VCF) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    
                    POS = int(fields[1])
                    REF = fields[3]
                    ALT =  fields[4].split(',')
                    
                    
                    # Get high-frequency allele if there are multiple ALT SNPs
                    if len(ALT) > 1:
                        info = fields[-1].split(':')
                        ad = [int(x) for x in info[1].split(',')]
                        ad_max = ad.index(max(ad))
                        ALT = ALT[ad_max - 1]
                    else:
                        ALT = ALT[0]
                        
                        
                    # Single nucleotide polymorphism
                    
                    if len(REF) == 1:
                        
                        # The ALT allele with the highest frequency is an insertion (quite rare ...): 
                        # skip, since first base of the insertion corresponds to the reference base
                        if len(ALT) > 1:
                            continue

                        if POS not in self.variants:
                            self.variants[POS] = {} 
                            self.variants[POS]['MTB_anc'] = REF  
                        
                        self.variants[POS][g] = ALT
                        

                    # Multi-nucleotide polymorphism
                    if len(REF) > 1 and len(REF) == len(ALT):

                        #print(g, fields, ALT)
                        self.stats['MNPs'] += 1

                        for i,base in enumerate(ALT):
                            if POS + i not in self.variants:
                                self.variants[POS+i] = {}
                                self.variants[POS+i]['MTB_anc'] = REF[i] 
                                
                            self.variants[POS+i][g] = ALT[i]
                            
 
    
        
        

class missing_positions:

    def __init__(self):
        
        
        self.depth_suffix = 'depth.gz'
        self.max_missing_site = 1
        self.max_missing_strain = 1
    
    def extract_rows(self, g_number):
        """
        Extracts specific rows from a file as a list of integers.
        
        :param file_path: Path to the file to read from.
        :param row_indices: List of row indices to extract.
        :return: List of missing positions for a given g_number
        """
        
        file_path = os.path.join(self.path_to_v2_output,g_number[0:3],g_number[3:5],g_number[5:],g_number+'.{}'.format(self.depth_suffix))
        variable_positions = list(self.variants.keys())
        
        with gzip.open(file_path, 'rb') as f:
            
            # Get the file size
            f.seek(0, 2)
            file_size = f.tell()
            
            # Initialize an empty list to store the extracted rows
            missing = []
            
            # Iterate over the row indices
            for row_idx in variable_positions:
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
                if row_int > self.minDepth:
                    missing.append(row_idx)
        
        return missing
    
    
    def process_chunk(self, chunk):
        """
        Processes a chunk of file_paths and row_indices in parallel, 
        using the extract_rows method.
        
        :param chunk: A list of tuples, where the first element of each tuple is a file path
            and the second element is a list of row indices to extract from that file.
        :return: A list of lists of missing positions. The outer list has the same length as the
            input chunk, and the inner lists have the same length as the number of row_indices
            in each element of the chunk.
        """
        results = []
        for g_number in chunk:
            results.append(self.extract_rows(g_number))
        return results
    
    
    def parallel_extract_rows(self, chunk_size = 100):
    
        """
        Extracts rows from a list of depth files in parallel.
        
        :param file_paths: A list of paths to depth files
        :param chunk_size: The number of depth files to process in parallel. Defaults to 1000.
        :return: A list of lists of missing positions. The outer list has the same length as the
            input list of depth files, and the inner lists have the same length as the number of
            row indices in each chunk.
        """
        chunks = [self.gnumbers[i:i + chunk_size] for i in range(0, len(self.gnumbers), chunk_size)]
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.nthreads) as executor:
            results = list(executor.map(self.process_chunk, chunks))
                       
        for g, missing in zip(self.gnumbers, results[0]):
            
            for pos in missing:
                                
                if g in self.variants[pos]:
                    self.sites_to_doublecheck['ALT allele and missing'] += (g, pos)
                
                else:
                    self.variants[pos][g] = '-'
                    self.sites_to_doublecheck['missing'] += (g, pos)
        
        
    def add_missing_bases(self, g_number, min_depth = 5):
        """ Traverse depth file for a given G number and record missing alleles at variable positions
        
        g_number : str
            G number
            
        variant_dictionary : dict
            Dictionary containing variant positions as keys.
            
        min_depth : int
            Positions with a coverage below this threshold will be considered missing

        """    
        
        missing = {}
    
        path_to_DEPTH =  os.path.join(self.path_to_v2_output,g_number[0:3],g_number[3:5],g_number[5:],g_number+'.{}'.format(self.depth_suffix))
        #path_to_DEPTH =  os.path.join(self.path_to_v2_output,g_number+'.{}'.format(self.depth_suffix))
        
        i = 1 # 1-based indexing in depth file
        
        #print('Traversing ' + g_number)
        
        with gzip.open(path_to_DEPTH,'rt') as f:

            for line in f:
                
                # Only consider missing bases at variable positions
                if i in self.variants:
                                   
                    fields = line.strip().split('\t')
                    depth = int(fields[0])
                    
                    if depth < min_depth:
    
                        #print('Low depth at ', i)
                        missing[i] = g_number
                        
                i += 1
                
        return missing
 

    def process_depth_files(self,  min_depth = 5):
               
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.nthreads) as executor:
            
            missing  = executor.map(self.add_missing_bases, self.gnumbers)
            
            for d in missing:
                for pos in d:
                    g = d[pos]
                    if g in self.variants[pos]:
                        self.sites_to_doublecheck['ALT allele and missing'] += (g, pos)
                 
                    else:
                        self.variants[pos][g] = '-'
            
            
    def add_outgroup_SNPs(self):
        """

        """
        pass



class write_output:

    def __init__(self, path_to_v2_output):
        pass

    def get_seqs(self, max_missing = 0.1):
        """ To do: evaluate site first, since at quite some sites only the
        reference differs.    
        
        """
        
        self.seqs = {g : '' for g in self.gnumbers}

        for pos in sorted(self.variants):
            
            site = ''
            
            if pos in self.repeats:
                self.stats['filt_repetitive'] += 1
                continue
            
            if pos in self.dr_loci:
                self.stats['filt_dr_loci'] += 1
                continue
                
            for g in self.gnumbers:
                
                # ALT ALLELE or missing
                if g in self.variants[pos]:
                    base = self.variants[pos][g]
                    
                    if len(base) != 1:
                        sys.stderr.write('More than 1 base!\n')
                        
                # REF ALLELE
                else:
                    base = self.variants[pos]['MTB_anc']
                        
                site += base


            # Evaluate site 
            if len(site) != len(self.gnumbers):
                self.sites_to_doublecheck['Number of alleles different from number of strains'] += (pos, )
                #sys.stdout.write(site + '\n')
                
            # Filter sites with too many missing bases
            n_missing = site.count('-')
            missing_prop = n_missing / len(self.gnumbers)
            if missing_prop > self.max_missing_site:
                self.stats['filt_site_missing'] += 1
                self.nr_sites_considered -= 1
                continue
                
            site_no_missing = site.replace('-', '')
            
            n_alleles = len(set(site_no_missing))
            
            if n_alleles == 1:
                #self.stats['only_ref_different'] += 1
                continue
            
            
            if n_alleles == 2:
                self.stats['biallelic'] += 1
            elif n_alleles > 2:
                self.stats['multiallelic'] += 1
                
            # add bases to sequences
            for base,g in zip(site, self.gnumbers):
                self.seqs[g] += base
                
            self.sites_in_alignment.append(pos)
                            
            
    def write_alignment(self, prefix, output_path):
        """ 
        
        """
        with open(output_path + '/' + prefix + '.aligned.fasta', 'w') as fasta_handle:
        
            for g in self.seqs:
                rec = SeqRecord(Seq(self.seqs[g]), id=g, name='', description='')
                SeqIO.write(rec, fasta_handle, 'fasta')

    
    def write_positions(self, prefix, output_path):
        with open(output_path + '/' + prefix + '.variable_sites.txt', 'w') as f:
              for pos in self.sites_in_alignment:
                  f.write(str(pos) + '\n')
            
    
    def write_stats(self, prefix, output_path):
        with open(output_path + '/' + prefix + '.alignment_stats.tsv', 'w') as f:
            for k in self.stats:
                row = k + '\t' + str(self.stats[k]) + '\n'
                f.write(row)
                
    def write_weird_sites(self):
        
        for k in self.sites_to_doublecheck:
            sys.stderr.write(k + '\n')
            for site in self.sites_to_doublecheck[k]:
                outline = '\t'.join(map(str,site) + '\n')
                sys.stderr.write(outline + '\n')
        
        pass
    
    