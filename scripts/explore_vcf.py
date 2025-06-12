#%%
import gzip

vcf = '/scicore/home/gagneux/stritt0001/TB/projects/MTBC_IS6110-poly/data/fixed_variants.merged.vcf.gz'

unique_genotypes = {}
n_alt_alleles = set()

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
        continue
    
    pos = int(fields[1])
    
    if pos > 696:
        break
    
    if pos == 696:
       
    
        # Get the REF and ALT alleles
        ref = fields[3]
        alt =  fields[4].split(',')
        alleles = [ref] + alt
        
        n_alt = len(alt)
        if n_alt not in n_alt_alleles:
            n_alt_alleles.add(len(alt))
        
        genotypes = fields[9:]
        for i, sample in enumerate(samples):
            fmt = genotypes[i].split(':')
            gt = fmt[0]
            if gt not in ['./.', '././.']:
                print(sample, gt, fmt, fields[:9])
        #%%
        if gt not in unique_genotypes:
            unique_genotypes[gt] = 0
            print(gt)
            
            if gt == '0/1/0':
                print(fields[:9], fmt, sample)
            
        unique_genotypes[gt] += 1