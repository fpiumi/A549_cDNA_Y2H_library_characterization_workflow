#!/bin/bash
module load bioinfo/NanoPlot/1.42.0

NanoPlot --fastq sequences.fastq.gz --tsv_stats --outdir nanoplot_output 
