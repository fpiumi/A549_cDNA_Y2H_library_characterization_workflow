#!/bin/bash
module load bioinfo/samtools/1.19

samtools sort -n -o tmp.qsort.bam tmp.bam