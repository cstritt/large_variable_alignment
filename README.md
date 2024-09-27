Create an alignment with variable positions in fasta format, using the 
output of the [TBEE variant calling pipeline](https://git.scicore.unibas.ch/TBRU/PipelineMTB/-/tree/master/Pipeline_TB_v2).

Problem when using many samples: traversing the file with the sequencing depth per positions takes ages with many samples  (20k+) ...

- Exclude repeats: https://git.scicore.unibas.ch/TBRU/PipelineMTB/-/tree/master/Pipeline_TB_v2/repetitive_regions
- Exclude DR loci: https://git.scicore.unibas.ch/TBRU/PipelineMTB/-/blob/master/Pipeline_TB/useful_scripts/20160911_DR_filter_pos_reseqTB.txt


# Optimize depth file traversal

