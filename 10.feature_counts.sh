#!/bin/bash
module load bioinfo/Subread/2.0.4

featureCounts -T 5 -t exon -g gene_id \
-a reference.gtf \
-o counts.txt minimap2_genome.sort.bam