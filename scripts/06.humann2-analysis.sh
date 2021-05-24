#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N humann2Analysis
#$ -pe shared 60
#$ -M obadbotanist@yahoo.com

# AUTHOR: Olabiyi Aderemi Obayomi
# EMAIL: obadbotanist@yahoo.com

# A script to perform humann2 analyis for function annotation using shogun/metatranscriptome reads while 
# simultaneously performing taxonomy profiling using metaphlan2

set -e

#activate the Metagenomics conda environment
source /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/bin/activate /gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics

export PERL5LIB='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/lib/site_perl/5.26.2/x86_64-linux-thread-multi' 

mpa_dir='/gpfs0/bioinfo/apps/Miniconda2/Miniconda_v4.3.21/envs/Metagenomics/bin'
PROJECT_DIR='/gpfs0/bioinfo/users/obayomi/metatranscriptomics/'
RAW_DATA_DIR="${PROJECT_DIR}/01.raw_data/"
OUT_DIR="${PROJECT_DIR}/11.Humann2_classify/"
UNIREF_DB='/gpfs0/bioinfo/databases/Humann2/uniref90/uniref'
CHOCOPHLAN_DB='/gpfs0/bioinfo/databases/Humann2/chocophlan/chocophlan'
SAMPLES=($(cat ${PROJECT_DIR}/metadata.tsv | awk 'NR>1{print $1}'))
# If no particular suffix was us just assign the empty string ''
SUFFIX=''
GENE_CLUSTER=90
INPUT_FORMAT='fastq.gz'
THREADS=10
# Path to metaplan_to_stamp.pl script of microbiome helper
METAPHLAN2STAMP='/gpfs0/bioinfo/users/obayomi/microbiome_helper/metaphlan_to_stamp.pl'

[ -d ${OUT_DIR} ] || mkdir ${OUT_DIR}


##on raw reads
parallel --jobs 5  "humann2 \
			--gap-fill on \
			--input-format ${INPUT_FORMAT} \
			--minpath on \
			--nucleotide-database ${CHOCOPHLAN_DB} \
			--protein-database ${UNIREF_DB} \
			--threads ${THREADS} \
			--output ${OUT_DIR}/{} \
			--input ${RAW_DATA_DIR}/{}/{}${SUFFIX}.${INPUT_FORMAT}"  ::: ${SAMPLES[*]}

# Adding code for normalizing genefamilies and pathabundance tables
parallel --jobs 5  "humann2_renorm_table \
			-i ${OUT_DIR}/{}/{}${SUFFIX}_genefamilies.tsv \
			-o ${OUT_DIR}/{}/{}${SUFFIX}_genefamilies.cpm.tsv \
			--units cpm"  ::: ${SAMPLES[*]}

parallel --jobs 5 "humann2_renorm_table \
			-i ${OUT_DIR}/{}/{}${SUFFIX}_pathabundance.tsv \
			-o ${OUT_DIR}/{}/{}${SUFFIX}_pathabundance.cpm.tsv \
			--units cpm" ::: ${SAMPLES[*]}


[ -d ${OUT_DIR}/bugs_list/ ] ||  mkdir ${OUT_DIR}/bugs_list/
[ -d ${OUT_DIR}/pathway_abund_tables/ ] ||  mkdir  ${OUT_DIR}/pathway_abund_tables/
[ -d  ${OUT_DIR}/gene_family_tables/ ] ||  mkdir  ${OUT_DIR}/gene_family_tables/
[ -d  ${OUT_DIR}/functional_groups_tables/ ] ||  mkdir  ${OUT_DIR}/functional_groups_tables/
[ -d  ${OUT_DIR}/renamed_groups_tables/ ] ||  mkdir  ${OUT_DIR}/renamed_groups_tables/


# copy the outfiles to one single folder for all the samples
function copy_files(){ 
	
	local SAMPLE=$1	
	local OUT_DIR=$2

	# metaphlan taxonomy tables
	cp  ${OUT_DIR}/${SAMPLE}/${SAMPLE}${SUFFIX}_humann2_temp/${SAMPLE}${SUFFIX}_metaphlan_bugs_list.tsv  ${OUT_DIR}/bugs_list/
	
	# pathway abundance tables
	cp  ${OUT_DIR}/${SAMPLE}/${SAMPLE}${SUFFIX}_pathabundance.tsv  ${OUT_DIR}/pathway_abund_tables/
	
	#gene families table
	cp  ${OUT_DIR}/${SAMPLE}/${SAMPLE}${SUFFIX}_genefamilies.tsv ${OUT_DIR}/gene_family_tables/

}


export -f copy_files
parallel --jobs 0 "copy_files {} ${OUT_DIR}" ::: ${SAMPLES[*]}

[ -d  ${OUT_DIR}/Unexported/raw/ ] || mkdir -p ${OUT_DIR}/Unexported/raw/
[ -d  ${OUT_DIR}/Unexported/relative_abundance/ ] || mkdir -p ${OUT_DIR}/Unexported/relative_abundance/
[ -d  ${OUT_DIR}/Unexported/metaphlan2/ ] || mkdir -p ${OUT_DIR}/Unexported/metaphlan2/

[ -d  ${OUT_DIR}/Export/raw/ ] || mkdir -p ${OUT_DIR}/Export/raw/
[ -d  ${OUT_DIR}/Export/relative_abundance/ ] || mkdir -p ${OUT_DIR}/Export/relative_abundance/
[ -d  ${OUT_DIR}/Export/metaphlan2/ ] || mkdir -p ${OUT_DIR}/Export/metaphlan2/


# Taxonomy profiling
# Merge metaphlan files
merge_metaphlan_tables.py ${OUT_DIR}/bugs_list/*  > ${OUT_DIR}/Unexported/metaphlan2/metaphlan2_merged.tsv

# To read this table into STAMP or R convert it be a .SPF file with this command:
${METAPHLAN2STAMP}  ${OUT_DIR}/Unexported/metaphlan2/metaphlan2_merged.tsv  > ${OUT_DIR}/Export/metaphlan2/metaphlan2_merged.spf
# Finally, notice that the sample names in our metadata don't match the column names in metaphlan2_merged.spf.
# These need to match for STAMP to read the table. This can be fixed by removing
# all instances of "_metaphlan_bugs_list" using the sed command.
sed -i "s/${SUFFIX}_metaphlan_bugs_list//g" ${OUT_DIR}/Export/metaphlan2/metaphlan2_merged.spf


# Functional profiling
# Pathways
# Raw counts (RPK)
# Join all individual humann2 tables tables into a single table
humann2_join_tables --input ${OUT_DIR}/pathway_abund_tables/ --file_name pathabundance --output  ${OUT_DIR}/Unexported/raw/humann2_pathabundance.tsv
humann2_split_stratified_table --input ${OUT_DIR}/Unexported/raw/humann2_pathabundance.tsv --output ${OUT_DIR}/Unexported/raw/ 
sed "s/${SUFFIX}_Abundance//g" ${OUT_DIR}/Unexported/raw/humann2_pathabundance_unstratified.tsv  >  ${OUT_DIR}/Export/raw/humann2_pathabundance_unstratified.spf
sed -i 's/# Pathway/MetaCyc_pathway/' ${OUT_DIR}/Export/raw/humann2_pathabundance_unstratified.spf

# Normalized counts
# Normalize each sample into relative abundance (so that the counts for each sample sum to 100)
humann2_renorm_table --input ${OUT_DIR}/Unexported/raw/humann2_pathabundance.tsv --units relab --output ${OUT_DIR}/Unexported/relative_abundance/humann2_pathabundance_relab.tsv
# Unstratify the pathway
humann2_split_stratified_table --input ${OUT_DIR}/Unexported/relative_abundance/humann2_pathabundance_relab.tsv --output ${OUT_DIR}/Unexported/relative_abundance/
sed "s/${SUFFIX}_Abundance//g" ${OUT_DIR}/Unexported/relative_abundance/humann2_pathabundance_relab_unstratified.tsv  >  ${OUT_DIR}/Export/relative_abundance/humann2_pathabundance_relab_unstratified.spf
sed -i 's/# Pathway/MetaCyc_pathway/' ${OUT_DIR}/Export/relative_abundance/humann2_pathabundance_relab_unstratified.spf



# Gene families
function regroup_tables(){	

	# function to regroup and rename gene families generated by humann2
	local SAMPLE=$1
	local OUT_DIR=$2
	local CLUSTER=$3
	local FUNCTIONAL_GROUPS=("uniref${CLUSTER}_go" "uniref${CLUSTER}_ko" \
                   "uniref${CLUSTER}_eggnog" "uniref${CLUSTER}_pfam" \
                   "uniref${CLUSTER}_level4ec" "uniref${CLUSTER}_infogo1000" \
                   "uniref${CLUSTER}_rxn")
	
	for group in ${FUNCTIONAL_GROUPS[*]}; do
			   
		group=$( echo ${group} | sed -E "s/uniref${CLUSTER}_//g")

		humann2_regroup_table \
			--input ${OUT_DIR}/${SAMPLE}/${SAMPLE}${SUFFIX}_genefamilies.tsv \
			--groups uniref${CLUSTER}_${group} \
			--output ${OUT_DIR}/functional_groups_tables/${SAMPLE}${SUFFIX}.${group}.tsv
	done

}


function add_names_to_tables(){

        # function to add descriptive names to regrouped gene family tabsles
        local SAMPLE=$1
        local OUT_DIR=$2
	local CLUSTER=$3
        local GROUPS_ONE=('infogo1000' 'go' 'pfam' 'eggnog') # names carried over from the regrouping step
	#local KEGG=('kegg-module' 'kegg-pathway' 'kegg-orthology')

        for group in ${GROUPS_ONE[*]}; do

                # Rename tables
                humann2_rename_table \
                        --input  ${OUT_DIR}/functional_groups_tables/${SAMPLE}${SUFFIX}.${group}.tsv \
                        --names ${group} \
                        --output ${OUT_DIR}/renamed_groups_tables/${SAMPLE}${SUFFIX}.${group}_names.tsv

        done

	# KeGG
	#for group in ${KEGG[*]}; do
	# Rename tables
        humann2_rename_table \
                  	--input  ${OUT_DIR}/functional_groups_tables/${SAMPLE}${SUFFIX}.ko.tsv \
                        --names kegg-orthology \
                        --output ${OUT_DIR}/renamed_groups_tables/${SAMPLE}${SUFFIX}.kegg-orthology_names.tsv
	#done

	# EC
	 humann2_rename_table \
                        --input  ${OUT_DIR}/functional_groups_tables/${SAMPLE}${SUFFIX}.level4ec.tsv \
                        --names ec \
                        --output ${OUT_DIR}/renamed_groups_tables/${SAMPLE}${SUFFIX}.ec_names.tsv

	# RXN
	humann2_rename_table \
                        --input  ${OUT_DIR}/functional_groups_tables/${SAMPLE}${SUFFIX}.rxn.tsv \
                        --names metacyc-rxn \
                        --output ${OUT_DIR}/renamed_groups_tables/${SAMPLE}${SUFFIX}.metacyc-rxn_names.tsv

	# uniref90
	humann2_rename_table \
                        --input  ${OUT_DIR}/${SAMPLE}/${SAMPLE}${SUFFIX}_genefamilies.tsv \
                        --names uniref${CLUSTER} \
                        --output ${OUT_DIR}/renamed_groups_tables/${SAMPLE}${SUFFIX}.uniref${CLUSTER}_names.tsv

}


export -f regroup_tables
parallel --jobs 0 "regroup_tables {}  ${OUT_DIR} ${GENE_CLUSTER}" ::: ${SAMPLES[*]}
export -f add_names_to_tables
parallel --jobs 0 "add_names_to_tables {}  ${OUT_DIR} ${GENE_CLUSTER}" ::: ${SAMPLES[*]}

function join_tables(){

        # function to regroup and rame gene families generated by humann2
        local OUT_DIR=$1
        local CLUSTER=$2
        local FUNCTIONAL_GROUPS=("uniref${CLUSTER}_go" "uniref${CLUSTER}_ko" \
                   "uniref${CLUSTER}_eggnog" "uniref${CLUSTER}_pfam" \
                   "uniref${CLUSTER}_level4ec" "uniref${CLUSTER}_infogo1000" \
                   "uniref${CLUSTER}_rxn")

        [ -d  ${OUT_DIR}/Unexported/raw/ ] || mkdir -p ${OUT_DIR}/Unexported/raw/
	[ -d  ${OUT_DIR}/Unexported/relative_abundance/ ] || mkdir -p ${OUT_DIR}/Unexported/relative_abundance/

        [ -d  ${OUT_DIR}/Export/raw/ ] || mkdir -p ${OUT_DIR}/Export/raw/
	[ -d  ${OUT_DIR}/Export/relative_abundance/ ] || mkdir -p ${OUT_DIR}/Export/relative_abundance/


	
        for group in ${FUNCTIONAL_GROUPS[*]}; do

                group=$( echo ${group} | sed -E "s/uniref${CLUSTER}_//g")

		# Raw tables
                # Join all individual humann2 ${group} tables tables into a single table
                humann2_join_tables --input  ${OUT_DIR}/functional_groups_tables/ --file_name ${group} --output  ${OUT_DIR}/Unexported/raw/humann2_${group}.tsv
		humann2_split_stratified_table --input ${OUT_DIR}/Unexported/raw/humann2_${group}.tsv --output ${OUT_DIR}/Unexported/raw/
		sed -E "s/${SUFFIX}_Abundance\-RPKs//g"  ${OUT_DIR}/Unexported/raw/humann2_${group}_unstratified.tsv  \
                >  ${OUT_DIR}/Export/raw/humann2_${group}_unstratified.spf
                sed -i -E "s/# Gene Family/${group}/"  ${OUT_DIR}/Export/raw/humann2_${group}_unstratified.spf
				
		# Normalized tables
                # Normalize each sample into relative abundance (so that the counts for each sample sum to 100)
                humann2_renorm_table --input  ${OUT_DIR}/Unexported/raw/humann2_${group}.tsv --units relab --output  ${OUT_DIR}/Unexported/relative_abundance/humann2_${group}_relab.tsv
                humann2_split_stratified_table --input ${OUT_DIR}/Unexported/relative_abundance/humann2_${group}_relab.tsv --output ${OUT_DIR}/Unexported/relative_abundance/
                # Using sed to make string replacements we can also format the header-line for input into STAMP.
                sed -E "s/${SUFFIX}_Abundance\-RPKs//g"  ${OUT_DIR}/Unexported/relative_abundance/humann2_${group}_relab_unstratified.tsv  \
                >  ${OUT_DIR}/Export/relative_abundance/humann2_${group}_relab_unstratified.spf
                sed -i -E "s/# Gene Family/${group}/"  ${OUT_DIR}/Export/relative_abundance/humann2_${group}_relab_unstratified.spf

	done

}

function join_named_tables(){

	local OUT_DIR=$1
	#local RENAME_GROUPS=('infogo1000' 'metacyc-rxn' 'kegg-module' \
        #                     'ec' 'go'  'pfam' 'eggnog' \
        #                     'kegg-pathway' 'kegg-orthology')
	local RENAME_GROUPS=('infogo1000' 'metacyc-rxn' \
                             'ec' 'go'  'pfam' 'eggnog' \
                             'kegg-orthology')

        [ -d  ${OUT_DIR}/Unexported/raw/ ] || mkdir -p ${OUT_DIR}/Unexported/raw/
	[ -d  ${OUT_DIR}/Unexported/relative_abundance/ ] || mkdir -p ${OUT_DIR}/Unexported/relative_abundance/

        [ -d  ${OUT_DIR}/Export/raw/ ] || mkdir -p ${OUT_DIR}/Export/raw/
	[ -d  ${OUT_DIR}/Export/relative_abundance/ ] || mkdir -p ${OUT_DIR}/Export/relative_abundance/



	for group in ${RENAME_GROUPS[*]}; do

                # Join all individual humann2 rennamed ${group} tables tables into a single table
                humann2_join_tables --input  ${OUT_DIR}/renamed_groups_tables/ --file_name ${group}_names --output ${OUT_DIR}/Unexported/raw/humann2_${group}_names.tsv
		humann2_split_stratified_table --input ${OUT_DIR}/Unexported/raw/humann2_${group}_names.tsv --output ${OUT_DIR}/Unexported/raw/

		 # Using sed to make string replacements we can also format the header-line for input into STAMP.
                sed -E  "s/${SUFFIX}_Abundance\-RPKs//g"  ${OUT_DIR}/Unexported/raw/humann2_${group}_names_unstratified.tsv  \
                >  ${OUT_DIR}/Export/raw/humann2_${group}_names_unstratified.spf
                sed -i "s/# Gene Family/${group}/"  ${OUT_DIR}/Export/raw/humann2_${group}_names_unstratified.spf
		

                # Normalize each sample into relative abundance (so that the counts for each sample sum to 100)
                humann2_renorm_table --input  ${OUT_DIR}/Unexported/raw/humann2_${group}_names.tsv --units relab --output  ${OUT_DIR}/Unexported/relative_abundance/humann2_${group}_names_relab.tsv
                humann2_split_stratified_table --input ${OUT_DIR}/Unexported/relative_abundance/humann2_${group}_names_relab.tsv --output ${OUT_DIR}/Unexported/relative_abundance/
                # Using sed to make string replacements we can also format the header-line for input into STAMP.
                sed -E  "s/${SUFFIX}_Abundance\-RPKs//g"  ${OUT_DIR}/Unexported/relative_abundance/humann2_${group}_names_relab_unstratified.tsv  \
                >  ${OUT_DIR}/Export/relative_abundance/humann2_${group}_names_relab_unstratified.spf
                sed -i "s/# Gene Family/${group}/"  ${OUT_DIR}/Export/relative_abundance/humann2_${group}_names_relab_unstratified.spf

        done

}


# Join and unstratitify unnamed tables
join_tables ${OUT_DIR} ${GENE_CLUSTER}
# Join and unstratitify named tables
join_named_tables ${OUT_DIR}

# Raw tables
# Join all individual humann2 gene_families tables tables into a single table 
humann2_join_tables --input  ${OUT_DIR}/gene_family_tables/ --file_name genefamilies --output   ${OUT_DIR}/Unexported/raw/humann2_genefamilies.tsv
humann2_split_stratified_table --input ${OUT_DIR}/Unexported/raw/humann2_genefamilies.tsv --output ${OUT_DIR}/Unexported/raw/
# Using sed to make string replacements we can also format the header-line for input into STAMP.
sed "s/${SUFFIX}_Abundance\-RPKs//g"  ${OUT_DIR}/Unexported/raw/humann2_genefamilies_unstratified.tsv  >  ${OUT_DIR}/Export/raw/humann2_genefamilies_unstratified.spf
sed -i 's/# Gene Family/Uniref_gene_family/'  ${OUT_DIR}/Export/raw/humann2_genefamilies_unstratified.spf

# Normalized tables
# Normalize each sample into relative abundance (so that the counts for each sample sum to 100)
humann2_renorm_table --input  ${OUT_DIR}/Unexported/raw/humann2_genefamilies.tsv --units relab --output  ${OUT_DIR}/Unexported/relative_abundance/humann2_genefamilies_relab.tsv
# Humann2 yields functions startified by taxa
# You could use grep to slice out only those pathways/genes that are linked to the genus of interest specifically e.g. grep  Streptococcus humanan2_stratified_file.tsv.
#In addition, more focused analyses could be based off of prior knowledge. For instance, pathways related to sugar degradation or vitamin B12 synthesis
# In order to focus on certain taxa or gene or pathway we need to to unstratify the pathway and gene families tables
humann2_split_stratified_table --input ${OUT_DIR}/Unexported/relative_abundance/humann2_genefamilies_relab.tsv --output ${OUT_DIR}/Unexported/relative_abundance/
# Then we can parse out the header line matching the gene or pathway of interest
# in this eaxmple (LACTOSECAT-PWY: lactose and galactose degradation I)
#head -n 1 humann2_pathabundance_relab_unstratified.tsv >> humann2_pathabundance_relab_LACTOSECAT-PWY.tsv
#grep "LACTOSECAT-PWY" humann2_pathabundance_relab_unstratified.tsv >> humann2_pathabundance_relab_LACTOSECAT-PWY.tsv
# Using sed to make string replacements we can also format the header-line for input into STAMP.
sed "s/${SUFFIX}_Abundance\-RPKs//g"  ${OUT_DIR}/Unexported/relative_abundance/humann2_genefamilies_relab_unstratified.tsv \
 >  ${OUT_DIR}/Export/relative_abundance/humann2_genefamilies_relab_unstratified.spf
sed -i 's/# Gene Family/Uniref_gene_family/'  ${OUT_DIR}/Export/relative_abundance/humann2_genefamilies_relab_unstratified.spf


# Cleaning Up
#[ -d  ${OUT_DIR}/Export ] || mkdir ${OUT_DIR}/Export
#[ -d  ${OUT_DIR}/Merged_tables ] || mkdir ${OUT_DIR}/Merged_tables
#mv ${OUT_DIR}/*.spf ${OUT_DIR}/Export  # files for export to R and Stamp
#mv ${OUT_DIR}/*.tsv ${OUT_DIR}/Merged_tables

# Read file into STAMP as before with the mapping file.
