#%%

import os 

#g_numbers = '/scicore/home/gagneux/stritt0001/TB/projects/MTBC_IS6110-poly/metadata/10k/gnumbers.txt'
g_numbers = '../testing/1k_genomes.txt'

path_to_v2 ='/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/v2'
vcf_suffix = 'mutect2.filtered.homo.snps.vcf'
depth_suffix = 'depth.gz'

prefix = '1k_genomes'



#%% Create list of depth files
outvcf = open(f'../testing/{prefix}.vcffiles.txt', 'w')
outdepth = open(f'../testing/{prefix}.depthfiles.tsv', 'w')

with open(g_numbers) as f:
    for line in f:
        g = line.strip()
        path_to_vcf =  os.path.join(path_to_v2,g[0:3],g[3:5],g[5:],g+'.{}'.format(vcf_suffix))
        path_to_depth =  os.path.join(path_to_v2,g[0:3],g[3:5],g[5:],g+'.{}'.format(depth_suffix))

        outvcf.write(path_to_vcf + '\n')
        outdepth.write(f'{g}\t{path_to_depth}\n')

outvcf.close()
outdepth.close()
        
# %%
