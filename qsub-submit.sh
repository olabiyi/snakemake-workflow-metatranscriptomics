#!/usr/bin/env bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N submit-jobs
#$ -pe shared 1

set -e

# Activate the main conda environment
source activate non_model_RNA_Seq
export PERL5LIB=/gpfs0/bioinfo/users/obayomi/miniconda3/envs/non_model_RNA_Seq/lib/site_perl/5.26.2/


# Generate the rule graph on the commadline
# Rule graph
# snakemake -s snakefile --rulegraph | dot -Tpng > rulegraph.png
# Directed Acyclic Graph (DAG)
# snakemake -s snakefile --dag | dot -Tpng > dag.png

# Run snmakemake on the cluster
# --jobs 100 # submit a maximum 100 jobs
# --latency-wait 60 # wait for 60 seconds before declaring that a job has failed
snakemake \
        --keep-going \
        --restart-times 3 \
        --rerun-incomplete  \
	--cluster-config config/config.yaml \
	--cluster 'qsub -q bioinfo.q -S /bin/bash -cwd -V -N {rule}.{wildcards} -e logs/{rule}/{rule}.{wildcards}.e -o logs/{rule}/{rule}.{wildcards}.o -pe shared {threads}' \
	--jobs 10 \
	--latency-wait 60 

