#!/bin/bash
module load bioinfo/samtools/1.19

samtools view minimap2_genome_sept24_k12.sam \
-bh \
-t /work/user/fpiumi/Data/RNA-Seq/Human/Homo_sapiens.GRCh38.cdna.all.fa \
-F 2324 > test_minimap_transcriptome.filt.bam

samtools sort test_minimap_transcriptome.filt.bam \
-o test_minimap_transcriptome.filt.sort.bam
samtools index test_minimap_transcriptome.filt.sort.bam

samtools view test_minimap_transcriptome.filt.sort.bam \
| cut -f 3 | sort | uniq -c