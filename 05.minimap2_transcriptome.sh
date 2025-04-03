#!/bin/bash
module load bioinfo/Minimap/2-2.26

minimap2 -ax map-ont -N 100 \
reference.cdna.fa \ 
sequences.fastq.gz > minimap2_transcriptome.sam
