# A549_cDNA_Y2H_library_characterization_workflow
The workflow describes a method to characterize a cDNA Y2H library using Nanopore sequencing


Step 1: genome indexation

Script name: 1.samtools_index_genome.sh


#!/bin/bash

module load bioinfo/samtools/1.19

samtools view faidx reference.fa


output : reference.fa.fai

indexed genome is in the same folder as the reference genome

