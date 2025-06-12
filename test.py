#%%
import os
import time
import sys

from lib import variable_alignment

class args:
    def __init__(self):
        
        dirpath = os.path.dirname(__file__)

        self.input = 'testing/1k_genomes.vcffiles.txt'
        self.output_folder = 'testing'
        
        self.depth_files = 'testing/1k_genomes.depthfiles.tsv'
        self.exclude = ''
        self.outvcf = os.path.join(dirpath,'resources/outgroup/G742339.mutect2.filtered.allvariants.AF_09.vcf.gz')
        self.outdepth = os.path.join(dirpath,'resources/outgroup/G742339.depth.gz')
        self.reference = os.path.join(dirpath,"resources/reference/MTB_ancestor_reference.fasta")

        self.mindepth = 5
        self.maxmissing = 1
        self.threads = 10
        
argumenti = args()


# Mise en place
mep = variable_alignment.mep(argumenti)

#%%
# Get variable positions
sys.stdout.write("Getting variable positions\n")
sys.stdout.flush()
start_varpos = time.time()
variantmatrix = variable_alignment.variantmatrix()
variantmatrix.add_SNPs(mep)
sys.stdout.write(f'{len(variantmatrix.variable_positions)} variable positions\n')
sys.stdout.flush()
end_varpos = time.time()
sys.stdout.write("Added variable positions in %f seconds\n" % (end_varpos - start_varpos))
sys.stdout.flush()


#%% 
# Convert to numpy array
start_numpy = time.time()
variantmatrix.convert_to_array(mep)
end_numpy = time.time()
sys.stdout.write("Converted to numpy array in %f seconds\n" % (end_numpy - start_numpy))



#%%
# Add missing positions and filter sites
sys.stdout.write("Adding missing positions\n")
sys.stdout.flush()
start_missing = time.time()
variantmatrix.traverse_depth_files_parallel(mep)
end_missing = time.time()
sys.stdout.write("Added missing positions in %f seconds\n" % (end_missing - start_missing))
sys.stdout.flush()
variantmatrix.max_missing_filter(mep)

#%%
# Write output
variantmatrix.write_output(mep)


# %%
