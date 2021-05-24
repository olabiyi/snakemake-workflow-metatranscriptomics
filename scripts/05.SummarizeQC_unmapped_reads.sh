#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N summarizeQC
#$ -pe shared 72

set -e

# A Script to qulity check sam and fastq files
# Here itis applied to summary the results of mapping reads to a host 
# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com
 

# Set the variables appropriately
source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

PROJECT_DIR='/gpfs0/bioinfo/projects/Amit_Gross/03.Metagenomics_stormwater_biofiltration/'
SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
SAM_SUFFIX='_mapped_and_unmapped'
FASTQ_DIR="${PROJECT_DIR}/03.remove_host/fastq_files"
FASTQ_SUFFIX='_host_removed.fastq.gz'
SAM_DIR="${PROJECT_DIR}/03.remove_host/"
SUMM_DIR="${PROJECT_DIR}/04.SummarizeQC_unmapped_reads/"

function process_sam_and_QC_reads(){

	local SAMPLE=$1
	local SAM_SUFFIX=$2
	local SAM_DIR=$3
	local SUMM_DIR=$4
	local FASTQ_SUFFIX=$5
	local FASTQ_DIR=$6
	
	[ -d ${SUMM_DIR}/${SAMPLE} ] || mkdir ${SUMM_DIR}/${SAMPLE}
	###########
	# Running samtools sort
	#----------------
	samtools \
		sort \
		--threads 10 \
		-o ${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam \
	 	${SAM_DIR}/${SAMPLE}${SAM_SUFFIX}.sam

	###########
	# Running samtools flagstat
	#----------------
	samtools \
		flagstat \
		${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam  \
		> ${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam.flagstat.txt

	###########
	# Running samtools index
	#----------------
	samtools \
		index \
		${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam

	###########
	# Running samtools stats
	#----------------
	samtools \
		stats \
		--remove-dups \
		${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam  \
		> ${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam.stats.txt

	###########
	# Running samtools idxstats
	#----------------
	samtools \
		idxstats \
		${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam  \
		> ${SUMM_DIR}/${SAMPLE}/${SAMPLE}${SAM_SUFFIX}.sorted.bam.idxstats.txt
	
	# Run fastqc
	fastqc \
	--threads 10 \
	--outdir ${SUMM_DIR}/${SAMPLE}/ \
	${FASTQ_DIR}/${SAMPLE}${FASTQ_SUFFIX}

}

# Quality to Check in parallel
export -f process_sam_and_QC_reads
parallel --jobs 10 "process_sam_and_QC_reads {} ${SAM_SUFFIX} ${SAM_DIR} ${SUMM_DIR} ${FASTQ_SUFFIX} ${FASTQ_DIR}" ::: ${SAMPLES[*]}

multiqc \
	--interactive \
	-s -f \
	-o ${SUMM_DIR} \
	${SUMM_DIR}

# Clean
rm -rf ${SAM_DIR}/*.sam
