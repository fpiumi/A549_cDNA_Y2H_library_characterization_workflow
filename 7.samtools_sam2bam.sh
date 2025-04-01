#!/bin/bash
module load bioinfo/samtools/1.19

samtools view -bS minimap2_genome.sam > minimap2_genome.bam