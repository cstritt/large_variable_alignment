# Create a SNP alignment

Create an alignment with variable positions in fasta format, using the 
output of the [TBEE variant calling pipeline](https://git.scicore.unibas.ch/TBRU/PipelineMTB/-/tree/master/Pipeline_TB_v2).


## Requirements
Biopython, numpy, joblib
```
conda create -n lva -c conda-forge biopython numpy joblib
```


## Input: 
* list of G numbers without header (-i)
* path to output directory (-o)

### For large alignments (> 10,000 strains):
* number of randomly chosen sites to include (-n)

## Defaults are provided for:
* drug resistance loci to exclude (resources/20160911_DR_filter_pos_reseqTB.txt)
* Illumina blindspots to exclude (resources/regions_blindspots_modlin_farhat.bed)
* reference genome (resources/MTB_ancestor_reference.fasta)
* outgroup strain (G742339, the recently described canettii strain ET1291)


## Differences comapred to previous alignment scripts
* each allele is recorded once, which present excessive memory use with many strains
* random access with the Python struct library is used to speed up the traversal of the depth files, which are used to identifiy missing sites
* an option for randomly subsampling sites. I couldn't figure out another way to handle many strains. Either memory requirements exploded or things took forever when creating an alignment for few hundred thousand variable sites. 



File origins: 
https://git.scicore.unibas.ch/TBRU/PipelineMTB/-/tree/master/Pipeline_TB_v2/repetitive_regions
https://git.scicore.unibas.ch/TBRU/PipelineMTB/-/blob/master/Pipeline_TB/useful_scripts/20160911_DR_filter_pos_reseqTB.txt


