#!/usr/bin/env bash

# Edit metadata

# Add a column that will contain the new names to be given to the files
paste -d '\t' \
	<(grep -v "#" old_metadata.tsv)  \
	<(awk 'NR>6' sample_name_file_Bacillus.txt) | \
	awk 'BEGIN{OFS="\t"} 
	{gsub("/gpfs0/biores/users/baubin/Noya/Transcriptome/", "00.raw_dir/",$6);
	gsub(".{1}$", "", $6); $7="00.raw_dir/"$1".fastq.gz" ;print($1,$2,$6,$7)}' \
	> metadata.tsv

# Add column names
(printf "%s\t%s\t%s\t%s\n" sample time old new; cat metadata.tsv)\
 > tmp && \
 mv tmp  metadata.tsv
 
# rename files in order to make them easy to work with
# old names
old=($(awk 'NR>1{print $3}' metadata.tsv)) 
# new names
new=($(awk 'NR>1{print $4}' metadata.tsv))
# rename files
for i in ${!old[*]}; do mv ${old[$i]} ${new[$i]};done

# put sample sequence in sample directory
mv 00.raw_dir/ 01.raw_data/
SAMPLES=($(awk 'NR > 1{print $1}' metadata.tsv))
cd 01.raw_data/
for sample in ${SAMPLES[*]}; do mkdir ${sample} && mv ${sample}.fastq.gz ${sample}/; done
