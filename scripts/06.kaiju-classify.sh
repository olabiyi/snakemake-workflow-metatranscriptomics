#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N classify_kaiju
#$ -pe shared 70
#$ -M obadbotanist@yahoo.com

# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

set -e

source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics
export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

PROJECT_DIR='/gpfs0/bioinfo/users/obayomi/metatranscriptomics/'
METADATA="${PROJECT_DIR}/metadata.tsv"
SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
KAIJU_PATH="/gpfs0/bioinfo/databases/kaiju_databases/kaijudb_archaea_bacteria_viral_euk_nr_2018"
KAIJU_CLASSIFY="${PROJECT_DIR}/09.Kaiju_classify/"
FMI="${KAIJU_PATH}/kaiju_db_nr_euk.fmi"
NODES="${KAIJU_PATH}/nodes.dmp"
NAMES="${KAIJU_PATH}/names.dmp"
RAW_DATA_DIR="${PROJECT_DIR}/01.raw_data/"
TAXON_LEVELS=(phylum class order family genus species)
FASTQ_SUFFIX='.fastq.gz'
KAIJU2TABLE="/gpfs0/bioinfo/apps/kaiju/kaiju_master_20190404/bin/kaiju2table"

#parallel --jobs 5 " [ -d ${KAIJU_CLASSIFY}/{} ] || mkdir ${KAIJU_CLASSIFY}/{} && kaiju \
#    -f ${FMI} \
#    -t ${NODES} \
#    -z 10 \
#    -E 1e-05 \
#    -i ${RAW_DATA_DIR}/{}/{}${FASTQ_SUFFIX} \
#    -o ${KAIJU_CLASSIFY}/{}/{}.kaiju.out" ::: ${SAMPLES[*]} 
	
# kaiju2krona - # Creating text report for krona and convert to html
parallel --jobs 0 "kaiju2krona \
        -n ${NAMES} \
        -t ${NODES} \
	-i ${KAIJU_CLASSIFY}/{}/{}.kaiju.out \
	-o ${KAIJU_CLASSIFY}/{}/{}.kaiju_out_4krona.txt && \
	ktImportText -o ${KAIJU_CLASSIFY}/{}/{}.kaiju.out.4krona.html \
	 ${KAIJU_CLASSIFY}/{}/{}.kaiju_out_4krona.txt " ::: ${SAMPLES[*]} 
  
# Create krona html chart report for all samples combined
find ${KAIJU_CLASSIFY}/ -type f -name "*_4krona.txt" |sort  > ${KAIJU_CLASSIFY}/krona_files.txt
FILES=($(find ${KAIJU_CLASSIFY}/ -type f -name "*_4krona.txt"))
basename -a -s '.kaiju_out_4krona.txt' ${FILES[*]} | sort  > ${KAIJU_CLASSIFY}/sample_names.txt
KTEXT_FILES=($(paste -d',' "${KAIJU_CLASSIFY}/krona_files.txt" "${KAIJU_CLASSIFY}/sample_names.txt"))
ktImportText -o "${KAIJU_CLASSIFY}/kaiju_report.html" ${KTEXT_FILES[*]}

# Kaiju2table
function run_kaiju2table() {

	local SAMPLE=$1
	local KAIJU_CLASSIFY=$2
	local KAIJU_PATH=$3
	TAXON_LEVELS=(phylum class order family genus species)
	NODES="${KAIJU_PATH}/nodes.dmp"
	NAMES="${KAIJU_PATH}/names.dmp"
	

	for TAXON_LEVEL in ${TAXON_LEVELS[*]}; do
		#echo ${TAXON_LEVEL}
		${KAIJU2TABLE} \
            		-t ${NODES} \
            		-n ${NAMES} \
            		-p  \
            		-r ${TAXON_LEVEL} \
            		-o ${KAIJU_CLASSIFY}/${SAMPLE}/${SAMPLE}.kaiju.summary.${TAXON_LEVEL}.tsv \
             		${KAIJU_CLASSIFY}/${SAMPLE}/${SAMPLE}.kaiju.out

	done

}

export -f run_kaiju2table
parallel --jobs 0 "run_kaiju2table {} ${KAIJU_CLASSIFY} ${KAIJU_PATH}" ::: ${SAMPLES[*]}

kaiju_out_FILES=($(find ${KAIJU_CLASSIFY} -type f -name "*.kaiju.out")) 

#merge tables
parallel --keep-order --jobs 0 "${KAIJU2TABLE} \
            -t ${NODES} \
            -n ${NAMES} \
            -p  \
            -r {} \
            -o ${KAIJU_CLASSIFY}/merged.kaiju.summary.{}.tsv \
             ${kaiju_out_FILES[*]}"  :::  ${TAXON_LEVELS[*]}

MERGED_FILES=($(find ${KAIJU_CLASSIFY} -type f -name "merged.kaiju.summary*"))

# covert the file names to sample names
parallel --jobs 0 "sed -i -E 's/.+\/(.+).kaiju\.out/\1/g' {} && sed -i -E 's/file/sample/' {}" ::: ${MERGED_FILES[*]


