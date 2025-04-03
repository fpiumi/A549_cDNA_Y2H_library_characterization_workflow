#!/bin/bash
module load bioinfo/Minimap/2-2.26

minimap2 -ax map-ont -N 100 \
reference.fa \
sequences.fastq.gz > minimap2_genome.sam