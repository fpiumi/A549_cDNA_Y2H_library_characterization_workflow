# A549_cDNA_Y2H_library_characterization_workflow
The workflow describes a method to characterize a cDNA Y2H library using Nanopore sequencing


## Step 1: genome indexation

Script name: 1.samtools_index_genome.sh
```shell
#!/bin/bash
module load bioinfo/samtools/1.19
samtools view faidx reference.fa
```
output : reference.fa.fai

indexed genome is in the same folder as the reference genome

## Step 2: Read quality assement (Nanoplot) 
Script name: 2.nanoplot.sh

```shell
#!/bin/bash
module load bioinfo/NanoPlot/1.42.0
NanoPlot --fastq sequences.fastq.gz
```

## Step 3: genome mapping
Script name: 3.minimap2_genome.sh

```shell
#!/bin/bash
module load bioinfo/Minimap/2-2.26
minimap2 -ax map-ont -N 100 \
reference.fa \
sequences.fastq.gz > minimap2_genome.sam
```
The –a argument specifies that the alignments should be output in SAM format (Sequence Alignment/Map), a standard format for representing sequence alignments. 
The -x map-ont preset adjusts minimap2's internal parameters to handle the characteristics of ONT reads, which are long reads (tens to hundreds of thousands of bases long) and which may have a higher error rate compared to other sequencing technologies, such as Illumina. This preset ensures robust alignments despite these errors. 
The -N parameter specifies the maximum number of secondary alignments reported per read. In genome alignment, the purpose is to map reads to the entire genomic sequence, which often contains repetitive elements, duplications, and homologous regions. Using -N 100 ensures to capture all possible mappings in a complex reference, especially in repetitive regions.

Command line:
sbatch minimap2_genome.sh
output : minimap2_genome.sam


## Step 4: splice-aware genome mapping
Script name: 4.minimap2_genome_splice.sh

```shell
#!/bin/bash
module load bioinfo/Minimap/2-2.26
minimap2 -ax splice –N 10 \
reference.fa \
sequences.fastq.gz > minimap2_genome_splice.sam
```
For splice-aware alignments with cDNA reads, the goal is typically to align reads to their true genomic origin while considering exon-exon junctions. Using a -N 10 ensures to reduce secondary alignments to focus on relevant mappings and limit alignments to biologically meaningful regions. 


## Step 5: transcriptome mapping
Script name: 5.minimap2_transcriptome.sh

```shell
#!/bin/bash
module load bioinfo/Minimap/2-2.26
minimap2 -ax map-ont –N 100 \
reference.cdna.fa \ 
sequences.fastq.gz > minimap2_transcriptome.sam
```

## Step 6: transcriptome and the non-coding transcriptome concatenation mapping

Concatenation (cat command) of cdna and ncrna references-> reference.cdna.ncrna.fa
Remove duplicate sequences, otherwise an error is raised during the BAM compression

```shell
awk '/^>/{f=!d[$1];d[$1]=1}f' reference.cdna.ncrna.fa > reference.cdna.ncrna.wo.duplicatesID.fa
```

Script name: 6.minimap2_coding_noncoding_transcriptome.sh

```shell
minimap2 -ax map-ont -N 100 \
reference.cdna.ncrna.wo.duplicatesID.fa \
sequences.fastq.gz > minimap2_transcriptome_coding_noncoding.sam
```



