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

## Step 6: transcriptome and non-coding transcriptome concatenation mapping

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


## Step 7: sam to bam compression 
Script name: 7.samtools_sam2bam.sh
```shell
#!/bin/bash
module load bioinfo/samtools/1.19
samtools view -bS minimap2_genome.sam > minimap2_genome.bam
```
This step has to be done for all the sam files. 
The -bS option in samtools view is used to specify that the input was in SAM format (-S) and that the output should be in BAM format (-b).

## Step 8: bam indexing/sorting
Script name: 8.samtools_sort_index.sh
```shell
#!/bin/bash
module load bioinfo/samtools/1.19
samtools sort minimap2_genome.bam -o minimap2_genome.sort.bam
samtools index minimap2_genome.sort.bam
```
This step has to be done for all the bam files. 



## Step 9: samtools flagstat
Script name: 9.samtools_flagstat.sh
```shell
#!/bin/bash
module load bioinfo/samtools/1.19
samtools flagstat minimap2_genome_sort.bam
```

samtools flagstat output result for the minimap2_genome_splice.sort.bam :
```
2017793 + 0 in total (QC-passed reads + QC-failed reads)
1237849 + 0 primary
617625 + 0 secondary
162319 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1732658 + 0 mapped (85.87% : N/A)
952714 + 0 primary mapped (76.97% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
How to interpret each line?
In sequence alignment, **primary**, **secondary**, and **supplementary alignments** are terms used to describe how reads align to a reference genome or sequence. These terms help categorize the relationships between a read and its alignment positions, particularly in cases where a read aligns to multiple locations or spans structural variations.

2017793 + 0 in total (QC-passed reads + QC-failed reads) : Total reads in the file (including primary, secondary, supplementary, and QC-failed reads). The first number represents QC-passed reads.The second number (after +) represents QC-failed reads.


1237849 + 0 primary : Primary alignments.
The primary alignment is the best alignment for a read, as determined by the aligner based on a scoring algorithm (e.g., highest alignment score, least mismatches). For paired-end reads, the primary alignment must also respect the pairing information.
In the SAM file:
- A primary alignment is indicated by the absence of the `0x100` flag (`FLAG` field in SAM format).
- There is only one primary alignment per read.

617625 + 0 secondary : Secondary alignments.
A secondary alignment represents an alternative alignment for a read, where the read aligns to multiple places in the reference sequence but is not considered the "best" alignment.
In the SAM file:
- A secondary alignment is marked with the `0x100` flag.
- Secondary alignments are typically output if the aligner detects multimapping reads (reads aligning to multiple locations) or when there is ambiguity in alignment.
- Secondary alignments do not split a read. The entire read is aligned at a different location.
If this number is high, you may have repetitive sequences or multimapped reads.

162319 + 0 supplementary : Supplementary alignments.
A supplementary alignment is used when a read spans large structural variations (e.g., translocations, inversions) or when it maps in a split manner across multiple regions of the reference. These alignments represent fragments of a single read that align to different parts of the reference.
In the SAM file:
- A supplementary alignment is marked with the `0x800` flag.
- The aligner typically uses supplementary alignments to represent split reads, where part of a read maps to one region and another part maps to a different region.
- Supplementary alignments are always tied to a primary alignment, and they form a complete representation of the read's alignment.
Typically seen in chimeric or split reads (e.g., structural variants, long-read sequencing).

0 + 0 duplicates : Reads marked as PCR duplicates (detected using samtools markdup or Picard). If high, it may indicate PCR amplification bias.

0 + 0 primary duplicates : The "duplicates" line counts all duplicate reads, including primary, secondary, and supplementary alignments.
The "primary duplicates" line counts only the primary alignments that are marked as duplicates.
This is useful because only primary reads are typically used for downstream analyses (e.g., variant calling).

1732658 + 0 mapped (85.87% : N/A) : Mapped reads: Reads that successfully aligned to the reference genome.
A high percentage (>90%) is good. If low, check for reference mismatches or contamination.

0 + 0 paired in sequencing : Total paired-end reads (if applicable).
If using single-end sequencing, this value should be 0.

0 + 0 read1 / 
0 + 0 read2 : 
The number of first and second reads in paired-end sequencing.

0 + 0 properly paired (N/A : N/A) : 
Reads correctly mapped in pairs (both reads align correctly with the expected insert size and orientation).
If this number is low, it might indicate low mapping quality or structural variations.

0 + 0 with itself and mate mapped : 
Both paired reads are mapped to the same reference sequence.

0 + 0 singletons (N/A : N/A) : 
Reads where only one mate is mapped.
A high number (>10%) could indicate low sequencing quality or a poor reference genome.

0 + 0 with mate mapped to a different chr : 
Reads whose mates map to a different chromosome.
Could suggest structural variations or contaminations.

0 + 0 with mate mapped to a different chr (mapQ>=5) : 
A more stringent count of reads mapped to different chromosomes, with a mapping quality (MAPQ) of at least 5.
If this number is high, it might indicate translocations or contamination.


## Step 10: Reads summarization
Script 10.feature_counts.sh
```shell
#!/bin/bash
module load bioinfo/Subread/2.0.4
featureCounts -T 5 -t exon -g gene_id \
-a reference.gtf \
-o minimap2_genome.counts.txt \
minimap2_genome.sort.bam
```
also working with splice-aware genome mapping
```shell
#!/bin/bash
module load bioinfo/Subread/2.0.4
featureCounts -T 5 -t exon -g gene_id 
-a reference.gtf \
 -o minimap2_genome_splice.counts.txt \
minimap2_genome_splice.sort.bam
```
Options :
-T specifies the number of threads to be used.
-t Specify the feature type. For example, features often correspond to exons and meta-features to genes
-g specify the attribute type used to group features (eg. exons) into meta-features. A meta-feature is the aggregation of a set of features. Here, meta-features are genes. The featureCounts program uses the gene_id attribute available in the GTF format annotation (or the GeneID column in the SAF format annotation) to group features into meta-features, ie. features belonging to the same meta-feature have the same gene identifier.
-a is the genome annotation file (gtf file).
-o specifies the name of the output file, which includes the read counts
file.bam is an alignment file: in this file, the reads we want to count are aligned to the same genome as the annotation file.

Output :
counts.txt
For each meta-feature, the “Length” column gives the total length of genomic regions covered by features included in that meta-feature. the “Length” column typically contains the total number of non-overlapping bases in exons belonging to the same gene for each gene.
Strand : exons number of a gene



## Step 11: Reads summarization with transcriptome bam files
feature counts is not working with transcriptome bam files, the following script must be used :

Script 11.samtools_counts_transcriptome_files
```shell
#!/bin/bash
module load bioinfo/samtools/1.19
samtools view minimap2_transcriptome.sam \
-bh \
-t reference.cdna.fa \
-F 2324 > minimap2_transcriptome.filt.bam

samtools sort minimap2_transcriptome.filt.bam \
-o minimap2_transcriptome.filt.sort.bam

samtools index minimap2_transcriptome.filt.sort.bam

samtools view minimap2_transcriptome.filt.sort.bam \
| cut -f 3 | sort | uniq -c
```

## Step 12: Count files annotation

An human genome annotation file is preliminarily downloaded from the biomart Ensembl Genes 112 database, (GRCh38.p14 human genes version) with the following attributes: Gene.stable.ID, Transcript.stable.ID, Gene.description, Gene.name, Gene.type 
```R
# annotation file opening
biomart <- read.csv2("mart_export.txt",sep="\t") 
library(dplyr)
biomart2 <- biomart %>% distinct(Gene.stable.ID, .keep_all = TRUE)

# counts file opening
counts_genome_splice <- read.table( "minimap2_genome_splice.counts.txt" ,sep="\t", header = TRUE ,na.strings = c("Gene.name" , "Gene.type" ) )
counts_genome_splice_short <- counts_genome_splice %>%
  dplyr::select(-c("Chr"   , "Start" ,  "End" ,  "Strand" , "Length"))

# variable name reduction: counts_genome_splice = GCS
GCS_annotated <- merge(counts_genome_splice_short, biomartbis, by.x="Geneid", by.y="Gene.stable.ID", all.x=TRUE)

# keep counts > 0
GCSA_above_zero <- GCS_annotated[which(GCS_annotated$minimap2_genome_splice.sort.bam > 0),]

# remove empty gene.names
CGSAAZ_not_empty <- GCSA_above_zero[!(GCSA_above_zero$Gene.name==""),]

# remove gene.names = NA
CGSAAZNE_wo_na <- na.exclude(CGSASupZ_not_empty)

# summarize counts
CGSAAZNE_wo_na_grouped <- CGSASupZNV_wo_na %>%
  group_by(Gene.name,Gene.type) %>%
  summarize(sum_Counts = sum(minimap2_genome_splice_sept24.sort.bam), #1
            na.rm = TRUE)  %>%
  arrange(Gene.name)

# remove duplicated gene.names
CGSASupZNV_wo_na_grouped2 <- CGSASupZNV_wo_na_grouped %>%
  group_by(Gene.name) %>%
  summarize(sum_Counts = sum(sum_Counts), #1
            na.rm = TRUE)  %>%
  arrange(Gene.name)

## add again annotation because a part of it was lost during the previous steps
biomartbis3 <- biomart2[!(biomart2$Gene.name==""),]

biomartbis4 <- biomartbis3 %>%
  group_by(Gene.name) %>%
  arrange(Gene.name)
CGSASupZNV_wo_na_grouped3 <- merge(CGSASupZNV_wo_na_grouped2,
                                   biomartbis3 , by.x="Gene.name")

```

## Step 13: Functional Analysis

Utilisation de ClusterProfiler sous R
```R
# Gene list obtention to run Clusterprofiler
gene_list <- CGSASupZNV_wo_na_grouped3$Gene.name

# groupGO
library(ClusterProfiler)
ggo_BP <- groupGO(gene     = gene_list,
                  OrgDb    = org.Hs.eg.db,
                  keyType  = 'SYMBOL',
                  ont      = "BP",
                  level    = 2,
                  readable = TRUE)
dim(ggo_BP)

## GO terms
ego_BP <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

barplot(ego_BP)

## GO terms reduction
library(GOSemSim)
library(enrichplot)

bp2 <- clusterProfiler::simplify(ego_BP, cutoff=0.7, by="p.adjust", select_fun=min)
d <- GOSemSim::godata('org.Hs.eg.db', ont="BP")
bp3 <- pairwise_termsim(bp2, method = "Wang", semData = d)
treeplot(bp3)


# Reactome terms
# get entrez IDs
mart <- biomaRt::useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl") #, 
#                   mirror = "useast")
getbm_res<-getBM(attributes = c("hgnc_symbol","entrezgene_id"), filters = "hgnc_symbol", 
                 values = CGSASupZNV_wo_na_grouped3$Gene.name, mart = mart, uniqueRows = F)
getbm_res$entrezgene_id<-as.character(getbm_res$entrezgene_id)
CGSASupZNV_wo_na_grouped3$entrezid <- getbm_res$entrezgene[match(CGSASupZNV_wo_na_grouped3$Gene.name,getbm_res$hgnc_symbol)]

library(ReactomePA)
packageVersion("ReactomePA")
reactome_res <- enrichPathway(gene=CGSASupZNV_wo_na_grouped2$entrezid, pvalueCutoff = 0.05, readable=TRUE)
```
