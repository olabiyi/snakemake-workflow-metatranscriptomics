#!/usr/bin/env bash

source activate non_model_RNA_Seq

export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/non_model_RNA_Seq/lib/5.26.2'

PROJECT_DIR='/gpfs0/bioinfo/users/obayomi/metatranscriptomics/'
PAIRED=False
METADATA="${PROJECT_DIR}/metadata.tsv"
SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))


if [ ${PAIRED}]
awk '{ if (NR%4==1) { gsub(/\s.*/,""); print $0"/1" } else { print } }' \
	/gpfs0/bioinfo/users/obayomi/non_model_RNA-Seq/test/data//samtools/Filter_unmapped_reads/C02.PL_FD1/C02.PL_FD1.bwa.view.sort.bam.F.fastq 


awk '{ if (NR%4==1) { gsub(/\s.*/,""); print $0"/1" } else { print } }' \
	/gpfs0/bioinfo/users/obayomi/non_model_RNA-Seq/test/data//samtools/Filter_unmapped_reads/C02.PL_FD1/C02.PL_FD1.bwa.view.sort.bam.F.fastq \
	> /gpfs0/bioinfo/users/obayomi/non_model_RNA-Seq/test/data//add_trinity_tags/Add_trinity_tags/C02.PL_FD1.bwa.view.sort.bam.F.fastq.trin_tags.fq
