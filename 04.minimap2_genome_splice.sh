#!/bin/bash
module load bioinfo/Minimap/2-2.26

minimap2 -ax splice -N 10 \
reference.fa \
sequences.fastq.gz > minimap2_genome_splice.sam