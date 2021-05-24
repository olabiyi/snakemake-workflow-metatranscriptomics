#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N remove_rna
#$ -pe shared 40
#$ -q bioinfo.q
#$ -V
#$ -cwd


set -e

source activate /gpfs0/bioinfo/users/obayomi/miniconda3/
export PERL5LIB=/gpfs0/bioinfo/users/obayomi/miniconda3/lib/5.32.0

SORTMERNA_DIR="/gpfs0/bioinfo/users/obayomi/metatranscriptomics/samsa2/programs/sortmerna-2.1"
SORTMERNA="${SORTMERNA_DIR}/sortmerna"
INDEX_DB="${SORTMERNA_DIR}/indexdb_rna"

REF_DIR="$SORTMERNA_DIR/rRNA_databases"
INDEX_DIR="$SORTMERNA_DIR/index"
PROJECT_DIR="/gpfs0/bioinfo/users/obayomi/metatranscriptomics"
THREADS=10
RAW_DATA_DIR="${PROJECT_DIR}/01.raw_data"
DATABASE="${REF_DIR}/silva-bac-16s-id90.fasta,${INDEX_DIR}/silva-bac-16s-db:${REF_DIR}/silva-bac-23s-id98.fasta,${INDEX_DIR}/silva-bac-23s-db:${REF_DIR}/silva-arc-16s-id95.fasta,${INDEX_DIR}/silva-arc-16s-db:${REF_DIR}/silva-arc-23s-id98.fasta,${INDEX_DIR}/silva-arc-23s-db:${REF_DIR}/silva-euk-18s-id95.fasta,${INDEX_DIR}/silva-euk-18s-db:${REF_DIR}/silva-euk-28s-id98.fasta,${INDEX_DIR}/silva-euk-28s-db:${REF_DIR}/rfam-5.8s-database-id98.fasta,${INDEX_DIR}/rfam-5.8s-database-db:${REF_DIR}/rfam-5s-database-id98.fasta,${INDEX_DIR}/rfam-5s-database-db"
UNZIP_DIR="${PROJECT_DIR}/04.unzip_fastqz"
DEPLETED_DIR="${PROJECT_DIR}/05.deplete_rrna"
INDEX_DATABASE='False' # Set to 'True' if the ribosomal rRNA databases have not been indexed
UNZIP_FASTQ='False' #'True' 
METADATA="${PROJECT_DIR}/metadata.tsv"
SAMPLES=($(awk 'NR>1{print $1}' ${METADATA}))

# index the refence databases
function index_DB(){
       
        local INDEX_DB=$1
        local REF_DIR=$2
        local INDEX_DIR=$3        

	${INDEX_DB} \
		-v \
		--ref ${REF_DIR}/silva-bac-16s-id90.fasta,${INDEX_DIR}/silva-bac-16s-db:${REF_DIR}/silva-bac-23s-id98.fasta,${INDEX_DIR}/silva-bac-23s-db:${REF_DIR}/silva-arc-16s-id95.fasta,${INDEX_DIR}/silva-arc-16s-db:${REF_DIR}/silva-arc-23s-id98.fasta,${INDEX_DIR}/silva-arc-23s-db:${REF_DIR}/silva-euk-18s-id95.fasta,${INDEX_DIR}/silva-euk-18s-db:${REF_DIR}/silva-euk-28s-id98.fasta,${INDEX_DIR}/silva-euk-28s-db:${REF_DIR}/rfam-5.8s-database-id98.fasta,${INDEX_DIR}/rfam-5.8s-database-db:${REF_DIR}/rfam-5s-database-id98.fasta,${INDEX_DIR}/rfam-5s-database-db  

}


function unzip_fastq(){

	local SAMPLE=$1
        local RAW_DATA_DIR=$2
        local OUTPUT_DIR=$3

	[ -d ${OUTPUT_DIR}/${SAMPLE}/ ] || mkdir -p  ${OUTPUT_DIR}/${SAMPLE}/

	# Unzip the fastq file because sortmerna 2.1 only works with unzipped files
	[ -e ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.fastq ] || \
	 zcat ${RAW_DATA_DIR}/${SAMPLE}/${SAMPLE}.fastq.gz \
	   > ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.fastq

}


function remove_rRNA(){
       
	local SAMPLE=$1
        local SORTMERNA=$2
        local DATABASE=$3
        local THREADS=$4
        local INPUT_DIR=$5
        local OUTPUT_DIR=$6 
	

        [ -d ${OUTPUT_DIR}/${SAMPLE}/ ] || mkdir -p ${OUTPUT_DIR}/${SAMPLE}/
	# remove 16S rRNA
	$SORTMERNA -a ${THREADS} \
    	--ref ${DATABASE} \
    	--reads ${INPUT_DIR}/${SAMPLE}/${SAMPLE}.fastq \
    	--aligned ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.ribosomes \
    	--other ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.ribodepleted \
    	--fastx \
    	--log \
    	-v

}

# Index ribosomal RNA databases
 
if [ ${INDEX_DATABASE} == 'True' ]; then  

	index_DB ${INDEX_DB} ${REF_DIR} ${INDEX_DIR}

fi


# Unzip fastq.gz files

if [ ${UNZIP_FASTQ} == 'True' ]; then

	export -f unzip_fastq

	parallel --jobs 4 "unzip_fastq {} ${RAW_DATA_DIR} ${UNZIP_DIR}"  ::: ${SAMPLES[*]}

fi


# Remove / deplete ribosomal RNA

export -f remove_rRNA

parallel --jobs 4 \
	"remove_rRNA {} ${SORTMERNA} ${DATABASE} ${THREADS} ${UNZIP_DIR} ${DEPLETED_DIR}"  \
	::: ${SAMPLES[*]} 

