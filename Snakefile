from os import path,getcwd

configfile: "config/config.yaml"
# used for generated merged kaiju tables
TAXON_LEVELS=['phylum', 'class', 'order', 'family', 'genus', 'species']
ruleorder: Diamond_blastx > Megan_classify
localrules: all, make_logs_directories,Generate_count_matrix

RULES = ["QC_pre_trim", "SummarizeQC_pre_trim", "Trim_reads", "QC_post_trim", "SummarizeQC_post_trim",
         "Build_rRNA_index", "Remove_rRNA", "QC_unmapped_reads", "SummarizeQC_unmapped", "Kaiju_classify",
         "Kaiju2krona", "Kaiju_ktImportText", "Kaiju2table", "Kaiju2_merge_tables", "Diamond_blastx",
         "Megan_classify", "Humann2_classify", "Metaphlan_merge", "Metaphlan2krona", "Metaphlan_ktextImport",
         "Humann2_group_tables", "Humann2_name_tables", "Humann2_copy_genes_and_pathways_tables",
         "Humann2_normalize", "Humann2_join_tables", "Humann2_unstratify_raw_tables","Humann2_unstratify_normalized_tables",
         "Add_trinity_tag", "Assemble_transcriptome", "prepare_transcriptome_reference","QC_Assembly", 
         "Generate_Gene_Transcript_Map","Estimate_abundance","Generate_count_matrix","Select_representative_transcript",
         "Generate_Gene_Transcript_Map_For_Selected", "Identify_coding_region","Predict_ribosomal_rna",
         "Swiss_prot_blastx", "Identify_protein_families", "Swiss_prot_blastp", "Generate_annotation_table"]

TARGET_RULES=["QC_pre_trim", "Trim_reads", "Build_rRNA_index", "prepare_transcriptome_reference", "Generate_annotation_table"]


# setting up datbase string for sortmerna
bac_16s=config['databases']['sortmerna']['bac']['16s']['ref'] + ','  + config['databases']['sortmerna']['bac']['16s']['index']
bac_23s=config['databases']['sortmerna']['bac']['23s']['ref'] + ','  + config['databases']['sortmerna']['bac']['23s']['index']
arch_16s=config['databases']['sortmerna']['arch']['16s']['ref'] + ','  + config['databases']['sortmerna']['arch']['16s']['index']
arch_23s=config['databases']['sortmerna']['arch']['23s']['ref'] + ','  + config['databases']['sortmerna']['arch']['23s']['index']
euk_18s=config['databases']['sortmerna']['euk']['18s']['ref'] + ','  + config['databases']['sortmerna']['euk']['18s']['index']
euk_28s=config['databases']['sortmerna']['euk']['28s']['ref'] + ','  + config['databases']['sortmerna']['euk']['28s']['index']
rfam_5s=config['databases']['sortmerna']['rfam']['5s']['ref'] + ','  + config['databases']['sortmerna']['rfam']['5s']['index']
rfam_5_8s=config['databases']['sortmerna']['rfam']['5_8s']['ref'] + ','  + config['databases']['sortmerna']['rfam']['5_8s']['index']

rRNA_DB=bac_16s + ":" + bac_23s + ":" + arch_16s + ":" + arch_23s + ":" + euk_18s + ":" + euk_28s + ":" + rfam_5s + ":" + rfam_5_8s

rule all:
    input:
        expand("logs/{RULE}/",RULE=TARGET_RULES),
        expand("02.QC/pre_trim/{sample}/{sample}_fastqc.html", 
               sample=config['samples']),
        "02.QC/pre_trim/multiqc_report.html",
        "04.QC/post_trim/multiqc_report.html",
        multiext(config['databases']['sortmerna']['rfam']['5s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['rfam']['5_8s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['bac']['16s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['bac']['23s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['arch']['16s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['arch']['23s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['euk']['18s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        multiext(config['databases']['sortmerna']['euk']['28s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
        expand("06.remove_rRNA/{sample}/{sample}.fastq.gz",
                sample=config['samples']),
        "07.QC/unmapped_reads/multiqc_report.html",
        "08.Kaiju_classify/{project}.html".format(project=config['project_name']),
        expand("08.Kaiju_classify/merged.{taxon_level}.tsv",
               taxon_level=TAXON_LEVELS),
        expand("08.Kaiju_classify/{sample}/{sample}.phylum.tsv",
               sample=config['samples']),
        expand("09.Megan_classify/{sample}/{sample}.tkn",
                sample=config['samples']),
        "10.Humann2_classify/Exported/metaphlan2/metaphlan2_merged.spf",
        "10.Humann2_classify/{project}.html".format(project=config['project_name']),
        "10.Humann2_classify/Unexported/raw/humann2_genefamilies_unstratified.tsv",
        "10.Humann2_classify/Unexported/relative_abundance/humann2_genefamilies_relab_unstratified.tsv",
        "12.Assemble_transcriptome/{project}_trinity/{project}_trinity.Trinity.fasta".format(project=config['project_name']),
	"13.QC_Assembly/report.html",
	"16.Generate_gene_count_matrix/{}.gene.counts.matrix".format(config['project_name']),
	"24.Generate_annotation_table/{}.trino_anno_rep.tsv".format(config['project_name'])


# This rule will make rule specific log directories
# in order to easily store the standard input and stand error
# generated when submiting jobs to the cluster
rule make_logs_directories:
    output:
        directory("logs/QC_pre_trim/"),
        directory("logs/Trim_reads/"),
        directory("logs/Build_rRNA_index/"),
        directory("logs/Humann2_unstratify_normalized_tables"),
        directory("logs/Generate_annotation_table/"),
        directory("logs/prepare_transcriptome_reference/")
    threads: 1
    shell:
        """
         [ -d logs/ ] || mkdir -p logs/
         cd logs/
         for RULE in {RULES}; do
          [ -d ${{RULE}}/ ] || mkdir -p ${{RULE}}/
         done
        """

# --------- Quality check, trimming and contaminants removal  -------------------

rule QC_pre_trim:
    input:
        read="01.raw_data/{sample}/{sample}.fastq.gz",
        log_dirs=rules.make_logs_directories.output
    output:
        "02.QC/pre_trim/{sample}/{sample}_fastqc.html"
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=5,
        o="logs/QC_pre_trim/{sample}/{sample}.log.o",
        e="logs/QC_pre_trim/{sample}/{sample}.log.e"
#    log: "logs/QC_pre_trim/{sample}/{sample}.log"
    threads: 5
    shell:
        "{params.program} --outdir {params.out_dir}/ "
        "--threads {params.threads} {input.read}"

rule SummarizeQC_pre_trim:
    input:
        expand(["02.QC/pre_trim/{sample}/{sample}_fastqc.html"],
                 sample=config['samples'])
    output:
        "02.QC/pre_trim/multiqc_report.html"
    params:
        program=config['programs_path']['multiqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        o="logs/SummarizeQC_pre_trim/multiqc.log.o",
        e="logs/SummarizeQC_pre_trim/multiqc.log.e"
    threads: 1
    shell:
        "{params.program} --interactive -f {params.out_dir} -o {params.out_dir}"


adaptors=config['parameters']['trimmomatic']['adaptors']
min_length=config['parameters']['trimmomatic']['min_len']
rule Trim_reads:
    input:
        read="01.raw_data/{sample}/{sample}.fastq.gz",
        log_dirs=rules.make_logs_directories.output
    output:
        "03.trimmed/{sample}/{sample}.fastq.gz"
#    log:
#        "logs/trimmomatic/{sample}/{sample}.log"
    params:
        program=config['programs_path']['trimmomatic'],
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length),
        o="logs/Trim_reads/{sample}/{sample}.log.o",
        e="logs/Trim_reads/{sample}/{sample}.log.e"
    threads: 10
    resources:
        mem_mb=1024
    shell:
        "{params.program} SE "
        "-threads {threads} "
        "{input.read} "
        "{output} "
        "{params.trimmer} > {log} 2>&1 "


rule QC_post_trim:
    input:
        "03.trimmed/{sample}/{sample}.fastq.gz"
    output:
        "04.QC/post_trim/{sample}/{sample}_fastqc.html"
    params:
        program=config['programs_path']['fastqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=5,
        o="logs/QC_post_trim/{sample}/{sample}.log.o",
        e="logs/QC_post_trim/{sample}/{sample}.log.e"
#    log: "logs/QC_post_trim/{sample}/{sample}.log"
    threads: 5
    shell:
        "{params.program} --outdir {params.out_dir} "
        "--threads {params.threads} {input} "


rule SummarizeQC_post_trim:
    input:
        expand("04.QC/post_trim/{sample}/{sample}_fastqc.html",
                 sample=config['samples'])
    output:
        "04.QC/post_trim/multiqc_report.html"
    params:
        program=config['programs_path']['multiqc'],
        out_dir=lambda w, output: path.dirname(output[0]),
        o="logs/SummarizeQC_post_trim/multiqc.log.o",
        e="logs/SummarizeQC_post_trim/multiqc.log.e"
#   log: "logs/SummarizeQC_post_trim/multiqc.log"
    threads: 1
    shell:
        "{params.program} --interactive -f {params.out_dir} -o {params.out_dir}"


rule Build_rRNA_index:
    input:
        multiext(config['databases']['sortmerna']['ref_dir'],
                      "rfam-5s-database-id98.fasta", "rfam-5.8s-database-id98.fasta",
                      "silva-bac-16s-id90.fasta", "silva-bac-23s-id98.fasta",
                      "silva-arc-16s-id95.fasta", "silva-arc-23s-id98.fasta",
                      "silva-euk-18s-id95.fasta", "silva-euk-28s-id98.fasta"),
        log_dirs=rules.make_logs_directories.output
    output:
        multiext(
            config['databases']['sortmerna']['rfam']['5s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['rfam']['5_8s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['bac']['16s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['bac']['23s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['arch']['16s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['arch']['23s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['euk']['18s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats"),
         multiext(
            config['databases']['sortmerna']['euk']['28s']['index'],
            ".bursttrie_0.dat", ".kmer_0.dat", ".pos_0.dat", ".stats")
        
#    log: "logs/Build_rRNA_index/Build_rRNA_index.log"
    threads: 10
    params:
        program=config['programs_path']['index_db'],
        o="logs/Build_rRNA_index/Build_rRNA_index.log.o",
        e="logs/Build_rRNA_index/Build_rRNA_index.log.e"
    shell:
        "{params.program} --ref {rRNA_DB} "


rule Remove_rRNA:
    input:
        rules.Build_rRNA_index.output,
        fastq="03.trimmed/{sample}/{sample}.fastq.gz"
    output:
        "06.remove_rRNA/{sample}/{sample}.fastq.gz"
    params:
        program=config['programs_path']['sortmerna'],
        threads=8,
        base_dir=lambda w, output: path.dirname(output[0]),
        o="logs/Remove_rRNA/{sample}/{sample}.log.o",
        e="logs/Remove_rRNA/{sample}/{sample}.log.e"
#    log: "logs/Remove_rRNA/{sample}/{sample}.log"
    threads: 10
    shell:
        """
        # unzip the gzip file because sortmena only 
        # works with unzipped files
        zcat {input.fastq} > {params.base_dir}/temp.fastq
        
        {params.program}  \
         --ref {rRNA_DB} \
         --reads {params.base_dir}/temp.fastq \
         --aligned  {params.base_dir}/{wildcards.sample}.ribosomes \
         --other {params.base_dir}/{wildcards.sample} \
         --log -a {threads} \
         --fastx && \
        # zip fastq
        gzip {params.base_dir}/{wildcards.sample}.fastq && \
        # clean up
        rm -rf {params.base_dir}/temp.fastq \
          {params.base_dir}/{wildcards.sample}.ribosomes.fastq 
        """

rule QC_unmapped_reads:
    input:
        "06.remove_rRNA/{sample}/{sample}.fastq.gz"
    output:
        "07.QC/unmapped_reads/{sample}/{sample}_fastqc.html"
    params:
        program = config['programs_path']['fastqc'],
        out_dir= lambda w, output: path.dirname(output[0]),
        threads=5,
        o="logs/QC_unmapped_reads/{sample}/{sample}.log.o",
        e="logs/QC_unmapped_reads/{sample}/{sample}.log.e"
#    log: "logs/QC_unmapped_reads/{sample}/{sample}.log"
    threads: 5
    shell:
        "{params.program} --outdir {params.out_dir} "
        "--threads {threads} {input}"


rule SummarizeQC_unmapped:
    input:
        expand("07.QC/unmapped_reads/{sample}/{sample}_fastqc.html",
               sample=config['samples'])
    output:
        "07.QC/unmapped_reads/multiqc_report.html"
    threads: 1
#    log: "logs/SummarizeQC_unmapped_reads/multiqc.log"
    params:
        program = config['programs_path']['multiqc'],
        out_dir = lambda w, output: path.dirname(output[0]),
        o="logs/SummarizeQC_unmapped_reads/multiqc.log.o",
        e="logs/SummarizeQC_unmapped_reads/multiqc.log.e"
    shell:
        "{params.program} --interactive "
        "-f {params.out_dir} -o {params.out_dir}"
 

# --------- Read-based Taxonomy classification ------------

#     ------------ Kaiju --------------

rule Kaiju_classify:
    input:
        "06.remove_rRNA/{sample}/{sample}.fastq.gz"
    output:
        "08.Kaiju_classify/{sample}/{sample}.kaiju.out"
    params:
        program=config['programs_path']['kaiju']['kaiju'],
        fmi=config['databases']['kaiju']['fmi'],
        nodes=config['databases']['kaiju']['nodes'],
        e_value=config['parameters']['kaiju']['evalue'],
        conda_activate=config['conda']['metagenomics']['env'],
        threads=10,
        o="logs/Kaiju_classify/{sample}/{sample}.log.o",
        e="logs/Kaiju_classify/{sample}/{sample}.log.e"        
#    log: "logs/Kaiju_classify/{sample}/{sample}.log"
    threads: 10
    shell:
        """
        set +u
       {params.conda_activate}
        set -u
       {params.program} -f {params.fmi} -t {params.nodes} \
        -z {params.threads} -E {params.e_value} \
        -i  {input} -o {output}
        """

rule Kaiju2krona:
    input:
        "08.Kaiju_classify/{sample}/{sample}.kaiju.out"
    output:
        "08.Kaiju_classify/{sample}/{sample}.krona.txt"
    params:
        program=config['programs_path']['kaiju']['kaiju2krona'],
        names=config['databases']['kaiju']['names'],
        nodes=config['databases']['kaiju']['nodes'],
        o="logs/Kaiju2krona/{sample}/{sample}.log.o",
        e="logs/Kaiju2krona/{sample}/{sample}.log.e"
#    log: "logs/Kaiju2krona/{sample}/{sample}.log"
    threads: 1
    shell:
        "{params.program} -i {input} -o {output} "
        " -t {params.nodes} -n {params.names} "


rule Kaiju_ktImportText:
    input:
        files=expand("08.Kaiju_classify/{sample}/{sample}.krona.txt",
                     sample=config['samples'])
    output:
        "08.Kaiju_classify/{project}.html".format(project=config['project_name'])
    params:
        program=config['programs_path']['kt_import_text'],
        o="logs/Kaiju_ktImportText/krona.log.o",
        e="logs/Kaiju_ktImportText/krona.log.e"
    threads: 10
#    log: "logs/Kaiju_ktImportText/krona.log"
    run:
        files = ["{file},{sample}".format(file=file, sample=sample)
                  for file,sample in zip(input.files,config['samples'])]
        shell("{params.program} -o {output} {files}")


rule Kaiju2table:
    input:
        "08.Kaiju_classify/{sample}/{sample}.kaiju.out"
    output:
        "08.Kaiju_classify/{sample}/{sample}.phylum.tsv",
        "08.Kaiju_classify/{sample}/{sample}.class.tsv",
        "08.Kaiju_classify/{sample}/{sample}.order.tsv",
        "08.Kaiju_classify/{sample}/{sample}.family.tsv",
        "08.Kaiju_classify/{sample}/{sample}.genus.tsv",
        "08.Kaiju_classify/{sample}/{sample}.species.tsv"
    params:
        program=config['programs_path']['kaiju']['kaiju2table'],
        names=config['databases']['kaiju']['names'],
        nodes=config['databases']['kaiju']['nodes'],
        prefix="08.Kaiju_classify/{sample}/{sample}",
        o="logs/Kaiju2table/{sample}/{sample}.log.o",
        e="logs/Kaiju2table/{sample}/{sample}.log.e"
#    log: "logs/Kaiju2table/{sample}/{sample}.log"
    threads: 2
    run:
        for taxon_level in TAXON_LEVELS:
            shell("{params.program} -t {params.nodes} -n {params.names} "
                  "-p -r {taxon_level} -o {params.prefix}.{taxon_level}.tsv "
                  "{input}")


rule Kaiju2_merge_tables:
    input:
        expand("08.Kaiju_classify/{sample}/{sample}.kaiju.out",
               sample=config['samples'])
    output:
        "08.Kaiju_classify/merged.phylum.tsv",
        "08.Kaiju_classify/merged.class.tsv",
        "08.Kaiju_classify/merged.order.tsv",
        "08.Kaiju_classify/merged.family.tsv",
        "08.Kaiju_classify/merged.genus.tsv",
        "08.Kaiju_classify/merged.species.tsv"
    params:
        program=config['programs_path']['kaiju']['kaiju2table'],
        names=config['databases']['kaiju']['names'],
        nodes=config['databases']['kaiju']['nodes'],
        prefix="08.Kaiju_classify/merged",
        o="logs/Kaiju2_merge_tables/merge.log.o",
        e="logs/Kaiju2_merge_tables/merge.log.e"
#    log: "logs/Kaiju2_merge_tables/merge.log"
    threads: 10
    run:
        for taxon_level in TAXON_LEVELS:
            shell("{params.program} -t {params.nodes} -n {params.names} "
                  "-p -r {taxon_level} -o {params.prefix}.{taxon_level}.tsv  {input}")
            shell("sed -i -E 's/.+\/(.+).kaiju\.out/\\1/g' {params.prefix}.{taxon_level}.tsv  && "
                  "sed -i -E 's/file/sample/' {params.prefix}.{taxon_level}.tsv")



# ------------------- Megan pipeline ------------------------------------

rule Diamond_blastx:
    input:
        "06.remove_rRNA/{sample}/{sample}.fastq.gz"
    output:
        "09.Megan_classify/{sample}/{sample}.daa"
    params:
        database=config['databases']['diamond']['nr'],
        program=config['programs_path']['diamond'],
        threads=10,
        o="logs/Diamond_blastx/{sample}/{sample}.log.o",
        e="logs/Diamond_blastx/{sample}/{sample}.log.e"
#    log: "logs/Diamond_blastx/{sample}/{sample}.log"
    threads: 10
    #threads: workflow.cores * 0.75 # set maximum threads to 75% of the avaliable cores
    resources:
        # remember to set the --restart-times option to how many attempy you
        # would like sankemake to try example --restart-times 3 for 3 attempts.
        # this will automatically increament the amount of memmory required
        # by 500MB for every failed attempt
        mem_mb=lambda wildcards, attempt: attempt * 500
    shell:
        "{params.program} blastx "
        "--query {input} --daa {output} "
        "--db {params.database} --threads {params.threads} --unal 1"

rule Megan_classify:
    input:
        # This just says that this rule is dependent on Diamond_blastx
        rules.Diamond_blastx.output
    output:
        temp("09.Megan_classify/{sample}/{sample}.tkn")
    params:
        program=config['programs_path']['megan']['meganizer'],
        ACC2TAXA=config['databases']['megan']['acc2taxa'],
        ACC2SEED=config['databases']['megan']['acc2seed'],
        ACC2INTERPRO=config['databases']['megan']['acc2interpro'],
        ACC2EGGNOG=config['databases']['megan']['acc2eggnog'],
        out_dir=lambda w, output: path.dirname(output[0]),
        o="logs/Megan_classify/{sample}/{sample}.log.o",
        e="logs/Megan_classify/{sample}/{sample}.log.e"
    threads: 10
#    log: "logs/Megan_classify/{sample}/{sample}.log"
    shell:
        """
        {params.program} \
        --in {input} \
        --minScore 50 \
        --maxExpected 0.01 \
        --minPercentIdentity 0 \
        --topPercent 5 \
        --minSupportPercent 0.05 \
        --minSupport 0 \
        --lcaAlgorithm weighted \
        --lcaCoveragePercent 80 \
        --readAssignmentMode readCount \
        --acc2taxa {params.ACC2TAXA} \
        --acc2seed {params.ACC2SEED} \
        --acc2interpro2go {params.ACC2INTERPRO} && \
        touch {params.out_dir}/{wildcards.sample}.tkn
        """        

# ---------- Humann2 Funtional classification and taxonomy assignment with Metaphlan2 -----
rule Humann2_classify:
    input:
        "06.remove_rRNA/{sample}/{sample}.fastq.gz"
    output:
        bugs_list="10.Humann2_classify/{sample}/{sample}_humann2_temp/"
                  "{sample}_metaphlan_bugs_list.tsv",
        gene_families="10.Humann2_classify/{sample}/{sample}_genefamilies.tsv",
        pathabundance="10.Humann2_classify/{sample}/{sample}_pathabundance.tsv",
        pathcoverage="10.Humann2_classify/{sample}/{sample}_pathcoverage.tsv"
    params:
        mpa_dir=config['programs_path']['humann2']['mpa_dir'],
        UNIREF_DB=config['databases']['humann2']['uniref90'],
        CHOCOPHLAN_DB=config['databases']['humann2']['chocophlan'],
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['humann2'],
        out_dir=lambda w, output: path.dirname(output.gene_families),
        o="logs/Humann2_classify/{sample}/{sample}.log.o",
        e="logs/Humann2_classify/{sample}/{sample}.log.e"
    threads: 10
#    log: "logs/Humann2_classify/{sample}/{sample}.log"
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} \
         --gap-fill on \
         --input-format fastq.gz \
         --minpath on \
         --nucleotide-database {params.CHOCOPHLAN_DB} \
         --protein-database {params.UNIREF_DB} \
         --threads {threads} \
         --output {params.out_dir} \
         --input {input}
        """

rule Metaphlan_merge:
    input:
        expand("10.Humann2_classify/{sample}/{sample}_humann2_temp/"
               "{sample}_metaphlan_bugs_list.tsv",
                sample=config['samples'])
    output:
        tsv="10.Humann2_classify/Unexported/metaphlan2/metaphlan2_merged.tsv",
        # STAMP format required for further downstream analyses
        spf="10.Humann2_classify/Exported/metaphlan2/metaphlan2_merged.spf"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['metaphlan']['merge_metaphlan'],
        metaplan2stamp=config['programs_path']['metaphlan']['metaphlan2stamp'],
        PERL5LIB=config['conda']['metagenomics']['perl5lib'],
        unexported_dir=lambda w, output: path.dirname(output.tsv),
        exported_dir=lambda w, output: path.dirname(output.spf),
        o="logs/Metaphlan_merge/merge.log.o",
        e="logs/Metaphlan_merge/merge.log.e"
    threads: 5
#    log: "logs/Metaphlan_merge/merge.log"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        # Make ouput direcetories if they don't already exist
        [ -d  {params.unexported_dir}/ ] || mkdir -p {params.unexported_dir}/
        [ -d  {params.exported_dir}/ ] || mkdir -p   {params.exported_dir}/

        {params.program} {input}  > {output.tsv}
        {params.metaplan2stamp} {output.tsv}  > {output.spf} && \
        # remove all instances of _metaphlan_bugs_list
        sed -i "s/_metaphlan_bugs_list//g" {output.spf}
        """

rule Metaphlan2krona:
    input:
        "10.Humann2_classify/{sample}/{sample}_humann2_temp/"
        "{sample}_metaphlan_bugs_list.tsv"
    output:
        "10.Humann2_classify/{sample}/{sample}_humann2_temp/{sample}.krona.txt"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['metaphlan']['metaphlan2krona'],
        o="logs/Metaphlan2krona/{sample}/{sample}.log.o",
        e="logs/Metaphlan2krona/{sample}/{sample}.log.e"
    threads: 1
#    log: "logs/Metaphlan2krona/{sample}/{sample}.log"
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} -p {input} -k {output}
        """


rule Metaphlan_ktextImport:
    input:
        files=expand("10.Humann2_classify/{sample}/{sample}_humann2_temp/{sample}.krona.txt",
                      sample=config['samples'])
    output:
        "10.Humann2_classify/{project}.html".format(project=config['project_name'])
    params:
        program=config['programs_path']['kt_import_text'],
        o="logs/Metaphlan_ktextImport/krona.log.o",
        e="logs/Metaphlan_ktextImport/krona.log.e"
    threads: 1
#    log: "logs/Metaphlan_ktextImport/krona.log"
    run:
        files = ["{file},{sample}".format(file=file, sample=sample)
                  for file,sample in zip(input.files,config['samples'])]
        shell("{params.program} -o {output} {files}")




# -------------- Function groups -----------------------------------------#

# Group the gene familes table to functional groups tables per sample

FUNCTION_GROUPS=["go", "ko","eggnog", "pfam", "level4ec", "infogo1000", "rxn"]

rule Humann2_group_tables:
    input: "10.Humann2_classify/{sample}/{sample}_genefamilies.tsv"
    output: 
        temp(expand("10.Humann2_classify/unnamed_function_groups_tables/"
               "{{sample}}/{{sample}}_{group}.tsv", group=FUNCTION_GROUPS)),
        genefamilies="10.Humann2_classify/unnamed_function_groups_tables/"
               "{sample}/{sample}_genefamilies.tsv"
#    log: "logs/Humann2_group_tables/{sample}/{sample}.log"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['regroup_tables'],
        out_dir=lambda w, output: path.dirname(output[0]),
        uniref="uniref90",
        o="logs/Humann2_group_tables/{sample}/{sample}.log.o",
        e="logs/Humann2_group_tables/{sample}/{sample}.log.e"
    threads: 5
    shell:
        """
        set +u
        {params.conda_activate}
        set -u

	for group in {FUNCTION_GROUPS};do

		{params.program} \
			--input {input} \
			--groups {params.uniref}_${{group}} \
			--output {params.out_dir}/{wildcards.sample}_${{group}}.tsv
        
	done
 
        # copy the gene family table to the unamed functions group directory
        # so that the next step may work
        cp {input} {output.genefamilies}
        
        """



# --------------- Add names to the function tables per sample --------------------

RENAME_GROUPS={"infogo1000": "infogo1000", "go": "go", 
               "pfam": "pfam", "eggnog": "eggnog",
               "ko": "kegg-orthology", "level4ec": "ec",
               "rxn": "metacyc-rxn", "genefamilies": "uniref90"}

# Add names to the function groups tables per sample
# The patways table comes already with names while the gene families
# don't come with names hence this step is not necessary for gene families
# and pathways
rule Humann2_name_tables:
    input:
        expand("10.Humann2_classify/unnamed_function_groups_tables/"
               "{{sample}}/{{sample}}_{group}.tsv", group=RENAME_GROUPS.keys())
    output:
        expand("10.Humann2_classify/named_function_groups_tables/"
               "{{sample}}_{name}_named.tsv",  name=list(RENAME_GROUPS.values()))
#    log: "logs/Humann2_name_tables/{sample}/{sample}.log"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['rename_tables'],
        out_dir=lambda w, output: path.dirname(output[0]),
        in_dir=lambda w, input: path.dirname(input[0]),
        o="logs/Humann2_name_tables/{sample}/{sample}.log.o",
        e="logs/Humann2_name_tables/{sample}/{sample}.log.e"
    threads: 5
    run:
        shell("set +u; {params.conda_activate}; set -u")

        for group, name in RENAME_GROUPS.items():

            shell("""
                  {program} \
                  --input  {in_dir}/{sample}_{group}.tsv \
                  --names {name} \
                  --output {out_dir}/{sample}_{name}_named.tsv
                  """.format(program=params.program, in_dir=params.in_dir,
                             sample=wildcards.sample,
                             group=group, out_dir=params.out_dir, name=name)
                 )



# ----------- Copy all the pathways and gene families to
# new directories for ease of joining the tables ------------------------------

rule Humann2_copy_genes_and_pathways_tables:
    input:
        gene_families=expand("10.Humann2_classify/{sample}/{sample}_genefamilies.tsv",
                      sample=config['samples']),
        pathabundance=expand("10.Humann2_classify/{sample}/{sample}_pathabundance.tsv", 
                       sample=config['samples'])
    output:
        pathabundance=expand("10.Humann2_classify/pathway_abund_tables/{sample}_pathabundance.tsv",
                sample=config['samples']),
        gene_families=expand("10.Humann2_classify/gene_family_tables/{sample}_genefamilies.tsv",
                sample=config['samples'])
    threads: 1
#    log: "logs/Humann2_copy_genes_and_pathways_tables/copy_tables.log.o"
    params:
        pathways_dir=lambda w, output: path.dirname(output.pathabundance[0]),
        genes_dir=lambda w, output: path.dirname(output.gene_families[0]),
        o="logs/Humann2_copy_genes_and_pathways_tables/copy_tables.log.o",
        e="logs/Humann2_copy_genes_and_pathways_tables/copy_tables.log.e"
    shell:
        """
        # copy pathways
        cp {input.pathabundance} {params.pathways_dir}
        # copy gene families
        cp {input.gene_families} {params.genes_dir} 
        """



# Join Humman2 tables
rule Humann2_join_tables:
    input: 
        gene_families=expand(rules.Humann2_copy_genes_and_pathways_tables.output.gene_families,
                     sample=config['samples']),
        pathabundance=expand(rules.Humann2_copy_genes_and_pathways_tables.output.pathabundance,
                     sample=config['samples']),
        unnamed_function_groups=expand(rules.Humann2_group_tables.output,
                     sample=config['samples']),
        named_function_groups=expand(rules.Humann2_name_tables.output,
                    sample=config['samples'])
    output: 
        pathabundance="10.Humann2_classify/Unexported/raw/"
                "humann2_pathabundance.tsv",
        gene_families="10.Humann2_classify/Unexported/"
                "raw/humann2_genefamilies.tsv",
        unnamed_function_groups=expand("10.Humann2_classify/Unexported/raw/without_names/"
                 "humann2_{group}.tsv", group=FUNCTION_GROUPS),
        named_function_groups=expand("10.Humann2_classify/Unexported/raw/with_names/"
                 "humann2_{name}_named.tsv", name=list(RENAME_GROUPS.values()))
#    log: "logs/Humann2_join_tables/Humann2_join_tables.log"
    threads: 5
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['join_tables'],
        pathway_dir=lambda w, input: path.dirname(input.pathabundance[0]),
        genes_dir=lambda w, input: path.dirname(input.gene_families[0]),
        unnamed_functions_dir=lambda w, input: path.dirname(input.unnamed_function_groups[0]),
        named_functions_dir=lambda w, input: path.dirname(input.named_function_groups[0]),
        unnamed_groups=FUNCTION_GROUPS,
        named_groups=list(RENAME_GROUPS.values()),
        unnamed_out_dir=lambda w, output: path.dirname(output.unnamed_function_groups[0]),
        named_out_dir=lambda w, output: path.dirname(output.named_function_groups[0]),
        o="logs/Humann2_join_tables/Humann2_join_tables.log.o",
        e="logs/Humann2_join_tables/Humann2_join_tables.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        set -u

        # Pathways tables
        {params.program} \
          --input {params.pathway_dir} \
          --file_name pathabundance \
          --output  {output.pathabundance}
         
        # Because the path abundance tables are already names
        # copy the table to the named folder
        cp {output.pathabundance} {params.named_out_dir}/

        # Gene families tables
        {params.program} \
          --input  {params.genes_dir} \
          --file_name genefamilies \
          --output  {output.gene_families}

                
        # Unamed function tables 
        INDIR=$(dirname  {params.unnamed_functions_dir})
        for group in {params.unnamed_groups}; do
        {params.program} \
          --input  ${{INDIR}} \
          --file_name ${{group}} \
          --search-subdirectories \
          --output  {params.unnamed_out_dir}/humann2_${{group}}.tsv
        done

        # Named function tables

        for group in {params.named_groups}; do
        {params.program} \
          --input  {params.named_functions_dir} \
          --file_name ${{group}}_named \
          --output  {params.named_out_dir}/humann2_${{group}}_named.tsv
        done

       """


# ---------------- Convert RPK to relative abundance ---------------------------------------

rule Humann2_normalize:
    input:
        gene_families=rules.Humann2_join_tables.output.gene_families,
        pathabundance=rules.Humann2_join_tables.output.pathabundance,
        unnamed_function_groups=rules.Humann2_join_tables.output.unnamed_function_groups,
        named_function_groups=rules.Humann2_join_tables.output.named_function_groups
    output:
        gene_families="10.Humann2_classify/Unexported/relative_abundance/"
                      "humann2_genefamilies_relab.tsv",
        pathabundance="10.Humann2_classify/Unexported/relative_abundance/"
                      "humann2_pathabundance_relab.tsv",
        unnamed_function_groups=expand("10.Humann2_classify/Unexported/relative_abundance/"
                       "without_names/humann2_{group}_relab.tsv", group=FUNCTION_GROUPS),
        named_function_groups=expand("10.Humann2_classify/Unexported/relative_abundance/"
                       "with_names/humann2_{name}_named_relab.tsv", name=list(RENAME_GROUPS.values()))
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['renorm'],
        unnamed_functions_in_dir=lambda w, input: path.dirname(input.unnamed_function_groups[0]),
        named_functions_in_dir=lambda w, input: path.dirname(input.named_function_groups[0]),
        unnamed_groups=FUNCTION_GROUPS,
        named_groups=list(RENAME_GROUPS.values()),
        unnamed_out_dir=lambda w, output: path.dirname(output.unnamed_function_groups[0]),
        named_out_dir=lambda w, output: path.dirname(output.named_function_groups[0]),
        o="logs/Humann2_normalize/Humann2_normalize.log.o",
        e="logs/Humann2_normalize/Humann2_normalize.log.e"
    threads: 1
#    log: "logs/Humann2_normalize/Humann2_normalize.log"
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
      
        # Gene families

        {params.program} \
         -i {input.gene_families} \
         -o {output.gene_families} \
         --units relab

        # Pathways

        {params.program} \
         -i {input.pathabundance} \
         -o {output.pathabundance} \
         --units relab

        
        # Unamed function tables

        for group in {params.unnamed_groups}; do
        {params.program} \
          -i  {params.unnamed_functions_in_dir}/humann2_${{group}}.tsv \
          --units relab \
          -o  {params.unnamed_out_dir}/humann2_${{group}}_relab.tsv
        done


        # Named function tables

        for group in {params.named_groups}; do
        {params.program} \
          -i  {params.named_functions_in_dir}/humann2_${{group}}_named.tsv \
          --units relab \
          -o  {params.named_out_dir}/humann2_${{group}}_named_relab.tsv
        done
        """


# Unstratify tables i.e remove species name from function assignments
# and rename the header of the exported tables

# Raw tables
rule Humann2_unstratify_raw_tables:
    input: 
        gene_families=rules.Humann2_join_tables.output.gene_families,
        pathabundance=rules.Humann2_join_tables.output.pathabundance,
        unnamed_function_groups=rules.Humann2_join_tables.output.unnamed_function_groups,
        named_function_groups=rules.Humann2_join_tables.output.named_function_groups
    output:
        unexported_gene="10.Humann2_classify/Unexported/raw/"
                        "humann2_genefamilies_unstratified.tsv",
        exported_gene="10.Humann2_classify/Exported/raw/"
                      "humann2_genefamilies_unstratified.spf",
        unexported_pathway="10.Humann2_classify/Unexported/raw/"
                           "humann2_pathabundance_unstratified.tsv",
        exported_pathway="10.Humann2_classify/Exported/raw/"
                         "humann2_pathabundance_unstratified.spf",
        unexported_unnamed_function_groups=expand("10.Humann2_classify/Unexported/raw/without_names/"
                 "humann2_{group}_unstratified.tsv", group=FUNCTION_GROUPS),
        unexported_named_function_groups=expand("10.Humann2_classify/Unexported/raw/with_names/"
                 "humann2_{name}_named_unstratified.tsv", name=list(RENAME_GROUPS.values())),
        exported_unnamed_function_groups=expand("10.Humann2_classify/Exported/raw/without_names/"
                 "humann2_{group}_unstratified.spf", group=FUNCTION_GROUPS),
        exported_named_function_groups=expand("10.Humann2_classify/Exported/raw/with_names/"
                 "humann2_{name}_named_unstratified.spf", name=list(RENAME_GROUPS.values()))
    threads: 10
#    log: "logs/Humann2_unstratify_raw_tables"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['split_stratified_table'],
        unexport_dir=lambda w, output: path.dirname(output.unexported_gene),
        export_dir=lambda w, output: path.dirname(output.exported_gene),
        unnamed_groups=FUNCTION_GROUPS,
        named_groups=list(RENAME_GROUPS.values()),
        o="logs/Humann2_unstratify_raw_tables/Humann2_unstratify_raw_tables.log.o",
        e="logs/Humann2_unstratify_raw_tables/Humann2_unstratify_raw_tables.log.e"
    shell:   
        """
        set +u
        {params.conda_activate}
        set -u


        # Gene families table
        {params.program} \
          --input {input.gene_families} \
          --output {params.unexport_dir}

        # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis
        sed "s/_Abundance\-RPKs//g"  {output.unexported_gene}  >  {output.exported_gene}
        sed -i 's/# Gene Family/Uniref_gene_family/'  {output.exported_gene}


        # Pathways table

        {params.program} \
          --input {input.pathabundance} \
          --output {params.unexport_dir}

        # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis

        sed "s/_Abundance//g" {output.unexported_pathway} >  {output.exported_pathway}
        sed -i 's/# Pathway/MetaCyc_pathway/' {output.exported_pathway}


        # Unnamed tables

       for group in {params.unnamed_groups}; do

           {params.program} \
             --input {params.unexport_dir}/without_names/humann2_${{group}}.tsv \
             --output {params.unexport_dir}/without_names/

            # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis.

            sed -E  "s/_Abundance\-RPKs//g"  {params.unexport_dir}/without_names/humann2_${{group}}_unstratified.tsv  \
             >  {params.export_dir}/without_names/humann2_${{group}}_unstratified.spf
            sed -i "s/# Gene Family/${{group}}/" {params.export_dir}/without_names/humann2_${{group}}_unstratified.spf
        

       done


        # Named tables

       for group in {params.named_groups}; do

           {params.program} \
             --input {params.unexport_dir}/with_names/humann2_${{group}}_named.tsv \
             --output {params.unexport_dir}/with_names/

           # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis.
           sed -E  "s/_Abundance\-RPKs//g"  {params.unexport_dir}/with_names/humann2_${{group}}_named_unstratified.tsv  \
            >  {params.export_dir}/with_names/humann2_${{group}}_named_unstratified.spf
           sed -i "s/# Gene Family/${{group}}/" {params.export_dir}/with_names/humann2_${{group}}_named_unstratified.spf

       done


        """



# Normalized tables
rule Humann2_unstratify_normalized_tables:
    input: 
        gene_families=rules.Humann2_normalize.output.gene_families,
        pathabundance=rules.Humann2_normalize.output.pathabundance,
        unnamed_function_groups=rules.Humann2_normalize.output.unnamed_function_groups,
        named_function_groups=rules.Humann2_normalize.output.named_function_groups
    output:
        unexported_gene="10.Humann2_classify/Unexported/relative_abundance/"
                        "humann2_genefamilies_relab_unstratified.tsv",
        exported_gene="10.Humann2_classify/Exported/relative_abundance/"
                      "humann2_genefamilies_relab_unstratified.spf",
        unexported_pathway="10.Humann2_classify/Unexported/relative_abundance/"
                           "humann2_pathabundance_relab_unstratified.tsv",
        exported_pathway="10.Humann2_classify/Exported/relative_abundance/"
                         "humann2_pathabundance_relab_unstratified.spf",
        unexported_unnamed_function_groups=expand("10.Humann2_classify/Unexported/relative_abundance/without_names/"
                 "humann2_{group}_relab_unstratified.tsv", group=FUNCTION_GROUPS),
        unexported_named_function_groups=expand("10.Humann2_classify/Unexported/relative_abundance/with_names/"
                 "humann2_{name}_named_relab_unstratified.tsv", name=list(RENAME_GROUPS.values())),
        exported_unnamed_function_groups=expand("10.Humann2_classify/Exported/relative_abundance/without_names/"
                 "humann2_{group}_relab_unstratified.spf", group=FUNCTION_GROUPS),
        exported_named_function_groups=expand("10.Humann2_classify/Exported/relative_abundance/with_names/"
                 "humann2_{name}_named_relab_unstratified.spf", name=list(RENAME_GROUPS.values()))
    threads: 10
#    log: "logs/Humann2_unstratify_normalized_tables"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['split_stratified_table'],
        unexport_dir=lambda w, output: path.dirname(output.unexported_gene),
        export_dir=lambda w, output: path.dirname(output.exported_gene),
        unnamed_groups=FUNCTION_GROUPS,
        named_groups=list(RENAME_GROUPS.values()),
        o="logs/Humann2_unstratify_normalized_tables/Humann2_unstratify_normalized_tables.log.o",
        e="logs/Humann2_unstratify_normalized_tables/Humann2_unstratify_normalized_tables.log.e"
    shell:   
        """
        set +u
        {params.conda_activate}
        set -u


        # Gene families table
        {params.program} \
          --input {input.gene_families} \
          --output {params.unexport_dir}

        # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis
        sed "s/_Abundance\-RPKs//g"  {output.unexported_gene}  >  {output.exported_gene}
        sed -i 's/# Gene Family/Uniref_gene_family/'  {output.exported_gene}


        # Pathways table

        {params.program} \
          --input {input.pathabundance} \
          --output {params.unexport_dir}

        # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis

        sed "s/_Abundance//g" {output.unexported_pathway} >  {output.exported_pathway}
        sed -i 's/# Pathway/MetaCyc_pathway/' {output.exported_pathway}


        # Unnamed tables

       for group in {params.unnamed_groups}; do

           {params.program} \
             --input {params.unexport_dir}/without_names/humann2_${{group}}_relab.tsv \
             --output {params.unexport_dir}/without_names/

            # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis.

            sed -E  "s/_Abundance\-RPKs//g"  {params.unexport_dir}/without_names/humann2_${{group}}_relab_unstratified.tsv  \
             >  {params.export_dir}/without_names/humann2_${{group}}_relab_unstratified.spf
            sed -i "s/# Gene Family/${{group}}/" {params.export_dir}/without_names/humann2_${{group}}_relab_unstratified.spf
        

       done


        # Named tables

       for group in {params.named_groups}; do

           {params.program} \
             --input {params.unexport_dir}/with_names/humann2_${{group}}_named_relab.tsv \
             --output {params.unexport_dir}/with_names/

           # Using sed to format the header-line to be used as input to R or STAMP for downstream analysis.
           sed -E  "s/_Abundance\-RPKs//g"  {params.unexport_dir}/with_names/humann2_${{group}}_named_relab_unstratified.tsv  \
            >  {params.export_dir}/with_names/humann2_${{group}}_named_relab_unstratified.spf
           sed -i "s/# Gene Family/${{group}}/" {params.export_dir}/with_names/humann2_${{group}}_named_relab_unstratified.spf

       done


        """




# ----------------------------- Assembly based analysis -----------------------------------------#
# Transcriptome assembly iusing Trinity and annotation using Trinotate 

# Add tags required by trinity for assembly
rule Add_trinity_tag:
    input: "06.remove_rRNA/{sample}/{sample}.fastq.gz"
    output: "11.Add_trinity_tag/{sample}/{sample}.tag.fastq"
    threads: 1
#    log: "logs/Add_trinity_tag/{sample}/{sample}.log"
    params:
        o="logs/Add_trinity_tag/{sample}/{sample}.log.o",
        e="logs/Add_trinity_tag/{sample}/{sample}.log.e"
    shell:
        """"
         awk '{{ if (NR%4==1) {{ gsub(/\s.*/,""); print $0"/1" }} else {{ print }} }}' \
        <(zcat {input}) > {output}
        """

# Assemble transcripts
# Trinity must form part of the output folder name else trinity will fail 
# Single end
rule Assemble_transcriptome:
    input: expand("11.Add_trinity_tag/{sample}/{sample}.tag.fastq", sample=config['samples'])
    output: "12.Assemble_transcriptome/{project}_trinity/"
            "{project}_trinity.Trinity.fasta".format(project=config['project_name'])
    threads: 50
#    log:"logs/Assemble_transcriptome/Assemble_transcriptome.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['trinity'],
        cpu=config['parameters']['trinity']['cpu'],
        max_memory=config['parameters']['trinity']['max_memory'],
        min_kmer_cov=config['parameters']['trinity']['min_kmer_cov'],
        o="logs/Assemble_transcriptome/Assemble_transcriptome.log.o",
        e="logs/Assemble_transcriptome/Assemble_transcriptome.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        
        FILES=$(echo {input} |sed -E 's/ /,/g')
        {params.program} \
         --CPU {params.cpu} \
         --max_memory {params.max_memory} \
         --seqType fq \
         --min_kmer_cov {params.min_kmer_cov} \
         --full_cleanup \
         --output {output} \
         --single ${{FILES}} 
        """

# Quality check the transcriptome assembly with quast
rule QC_Assembly:
    input: rules.Assemble_transcriptome.output
    output: "13.QC_Assembly/report.html"
    threads: 10
#    log: "logs/QC_Assembly/QC_Assembly.log"
    params:
        program=config['programs_path']['quast'],
        threads=10,
        out_dir=lambda w, output: path.dirname(output[0]),
        o="logs/QC_Assembly/QC_Assembly.log.o",
        e="logs/QC_Assembly/QC_Assembly.log.e"
    shell:
        "{params.program} "
        "-t {params.threads} "
        "--output-dir {params.out_dir}  {input}"

# -------------------- Estimate transcript abundance -------------------------#

rule Generate_Gene_Transcript_Map:
    input: rules.Assemble_transcriptome.output
    output: "14.Generate_Gene_Transcript_Map/{}_trinity"
            ".Trinity.gene_trans_map".format(config['project_name'])
    threads: 10
#    log:"logs/Assemble_transcriptome/Assemble_transcriptome.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['gene_trans_map'],
        o="logs/Generate_Gene_Transcript_Map/Generate_Gene_Transcript_Map.log.o",
        e="logs/Generate_Gene_Transcript_Map/Generate_Gene_Transcript_Map.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} {input} > {output}
        """

# prepare bowtie reference for reads count estimation
rule prepare_transcriptome_reference:
    input:
        transcriptome=rules.Assemble_transcriptome.output,
        gene_trans_map=rules.Generate_Gene_Transcript_Map.output
    output:
        multiext(
        "12.Assemble_transcriptome/{project}_trinity/"
        "{project}_trinity.Trinity.fasta".format(project=config['project_name']),
        ".bowtie2.1.bt2", ".bowtie2.2.bt2", ".bowtie2.rev.1.bt2", ".bowtie2.rev.2.bt2")
    threads: 10
#    log: "logs/prepare_reference/prepare_reference.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['estimate_abundance'],
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=10,
        o="logs/prepare_reference/prepare_reference.log.o",
        e="logs/prepare_reference/prepare_reference.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
        --est_method RSEM \
        --prep_reference \
        --gene_trans_map  {input.gene_trans_map} \
        --thread_count {params.threads} \
        --aln_method bowtie2 \
        --transcripts {input.transcriptome}
        """

rule Estimate_abundance:
    input:
        rules.prepare_transcriptome_reference.output,
        read="11.Add_trinity_tag/{sample}/{sample}.tag.fastq",
        transcriptome=rules.Assemble_transcriptome.output,
        gene_trans_map=rules.Generate_Gene_Transcript_Map.output
    output:
        "15.Estimate_abundance/{sample}/RSEM.isoforms.results",
        "15.Estimate_abundance/{sample}/RSEM.genes.results"
    threads: 10
#    log:"logs/Estimate_abundance/{sample}.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['estimate_abundance'],
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=10,
        o="logs/Estimate_abundance/{sample}.log.o",
        e="logs/Estimate_abundance/{sample}.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
        --est_method RSEM \
        --gene_trans_map  {input.gene_trans_map} \
        --seqType fq \
        --thread_count {params.threads} \
        --aln_method bowtie2 \
        --transcripts {input.transcriptome} \
        --output_dir {params.out_dir} \
        --single {input.read}
        """

rule Generate_count_matrix:
    input: 
        isoforms=expand("15.Estimate_abundance/{sample}/RSEM.isoforms.results",
                        sample=config['samples']),
        gene_trans_map=rules.Generate_Gene_Transcript_Map.output
    output:
        gene_matrix="16.Generate_gene_count_matrix/{}"
                    ".gene.counts.matrix".format(config['project_name']),
        gene_tpm="16.Generate_gene_count_matrix/{}"
                    ".gene.TPM.not_cross_norm".format(config['project_name']),
        isoform_matrix="16.Generate_gene_count_matrix/{}"
                       ".isoform.counts.matrix".format(config['project_name']),
        isoform_tpm="16.Generate_gene_count_matrix/{}"
                    ".isoform.TPM.not_cross_norm".format(config['project_name'])
    threads: 10
#    log:"logs/Generate_count_matrix/Generate_count_matrix.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['count_matrix'],
        prefix="16.Generate_gene_count_matrix/{}".format(config['project_name']),
        threads=10,
        o="logs/Generate_count_matrix/Generate_count_matrix.log.o",
        e="logs/Generate_count_matrix/Generate_count_matrix.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
         --est_method RSEM \
         --name_sample_by_basedir \
         --gene_trans_map  {input.gene_trans_map} \
         --out_prefix {params.prefix} \
         {input.isoforms}
        """

# ----------- Annotation of selected (highly expressed) transcripts using Trinotate ----------#

rule Select_representative_transcript:
    input: 
        isoform_tpm=rules.Generate_count_matrix.output.isoform_tpm,
        transcriptome=rules.Assemble_transcriptome.output,
        gene_trans_map=rules.Generate_Gene_Transcript_Map.output
    output:
        "17.Select_representative_transcript/{}"
        ".assembly.fasta".format(config['project_name'])
    threads: 10
#    log:"logs/Generate_count_matrix/Generate_count_matrix.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['filter_transcript'],
        prefix="16.Generate_gene_count_matrix/{}".format(config['project_name']),
        threads=10,
        o="logs/Select_representative_transcript/Select_representative_transcript.log.o",
        e="logs/Select_representative_transcript/Select_representative_transcript.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
         --highest_iso_only \
         --matrix {input.isoform_tpm} \
         --transcripts {input.transcriptome} \
         --gene_to_trans_map {input.gene_trans_map} \
         > {output}
        """

rule Generate_Gene_Transcript_Map_For_Selected:
    input: rules.Select_representative_transcript.output
    output: "18.Generate_Gene_Transcript_Map_For_Selected/{}"
            ".assembly.fasta.gene_trans_map".format(config['project_name'])
    threads: 10
#    log:"logs/Generate_Gene_Transcript_Map_For_Selected/"
#        "Generate_Gene_Transcript_Map_For_Selected.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinity']['gene_trans_map'],
        o="logs/Generate_Gene_Transcript_Map_For_Selected/"
           "Generate_Gene_Transcript_Map_For_Selected.log.o",
        e="logs/Generate_Gene_Transcript_Map_For_Selected/"
          "Generate_Gene_Transcript_Map_For_Selected.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} {input} > {output}
        """

# Identify protein coding sequences using Transdecoder
rule Identify_coding_region:
    input: 
        transcriptome=rules.Select_representative_transcript.output,
        gene_trans_map=rules.Generate_Gene_Transcript_Map_For_Selected.output
    output: "19.Identify_coding_region/longest_orfs.pep"
    threads: 10
#    log:"logs/Identify_coding_region/Identify_coding_region.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinotate']['transdecoder'],
        out_dir=lambda w, output: path.dirname(output[0]),
        o="logs/Identify_coding_region/Identify_coding_region.log.o",
        e="logs/Identify_coding_region/Identify_coding_region.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
           --output_dir {params.out_dir} \
           -t {input.transcriptome} \
           --gene_trans_map {input.gene_trans_map}
        """

# Predict ribosomal RNA using RNAmmer - spewcify organism type as bacteria
rule Predict_ribosomal_rna:
    input: rules.Select_representative_transcript.output
    output: "20.Predict_ribosomal_rna/{}.assembly"
            ".fasta.rnammer.gff".format(config['project_name'])
    threads: 10
#    log:"logs/Predict_ribosomal_rna/Predict_ribosomal_rna.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinotate']['rnammer']['rnammer'],
        rnammer_path=config['programs_path']['trinotate']['rnammer']['path'],
        out_dir=lambda w, output: path.dirname(output[0]),
        trans_full_path=lambda w, input: getcwd() + "/" + input[0],
        o="logs/Predict_ribosomal_rna/Predict_ribosomal_rna.log.o",
        e="logs/Predict_ribosomal_rna/Predict_ribosomal_rna.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
        
        cd {params.out_dir}

        {params.program} \
          --path_to_rnammer {params.rnammer_path} \
          --transcriptome {params.trans_full_path} \
          --org_type bac 
        """

rule Swiss_prot_blastx:
    input: rules.Select_representative_transcript.output
    output: "21.Swiss_prot_blastx/{}.assembly"
            ".blast.out".format(config['project_name'])
    threads: 40
#    log:"logs/Swiss_prot_blastx/Swiss_prot_blastx.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinotate']['blastx'],
        parallel=config['programs_path']['parallel'],
        threads=40,
        database=config['databases']['trinotate']['sprot'],
        o="logs/Swiss_prot_blastx/Swiss_prot_blastx.log.o",
        e="logs/Swiss_prot_blastx/Swiss_prot_blastx.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
           -outfmt 6 \
           -max_target_seqs 1 \
           -num_threads {params.threads} \
           -db {params.database} \
           -query {input} > {output}
        """

# Predict PFam protein families using hmm scan from the predicted coding sequences
rule Identify_protein_families:
    input: rules.Identify_coding_region.output
    output: "22.Identify_protein_families/{}.assembly"
            ".hmmscan.domtblout".format(config['project_name'])
    threads: 30
    log:"logs/Identify_protein_families/Identify_protein_families.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinotate']['hmmscan'],
        threads=30,
        database=config['databases']['trinotate']['pfam'],
        o="logs/Identify_protein_families/Identify_protein_families.log.o",
        e="logs/Identify_protein_families/Identify_protein_families.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        {params.program} \
          --cpu {params.threads} \
          --domtblout {output} \
          {params.database} \
          {input} > {log} 2>&1          
         """ 

rule Swiss_prot_blastp:
    input: rules.Identify_coding_region.output
    output: "23.Swiss_prot_blastp/{}.assembly"
            ".blast.out".format(config['project_name'])
    threads: 40
#    log:"logs/Swiss_prot_blastp/Swiss_prot_blastp.log"
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinotate']['blastp'],
        parallel=config['programs_path']['parallel'],
        threads=2,
        database=config['databases']['trinotate']['sprot'],
        o="logs/Swiss_prot_blastp/Swiss_prot_blastp.log.o",
        e="logs/Swiss_prot_blastp/Swiss_prot_blastp.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u
    
        cat {input} | \
        {params.parallel} \
          --jobs 20 \
          --block 100k \
          --recstart '>' \
          --pipe {params.program} \
           -outfmt 6 \
           -max_target_seqs 1 \
           -num_threads {params.threads} \
           -db {params.database} \
           -query - > {output}
        """

rule Generate_annotation_table:
    input:
        transcriptome=rules.Select_representative_transcript.output,
        gene_trans_map=rules.Generate_Gene_Transcript_Map_For_Selected.output, 
        coding_sequences=rules.Identify_coding_region.output,
        swissprot_blastx=rules.Swiss_prot_blastx.output,
        swissprot_blastp=rules.Swiss_prot_blastp.output,
        pfam=rules.Identify_protein_families.output,
        rRNA=rules.Predict_ribosomal_rna.output
    output: 
        sqlitedb="24.Generate_annotation_table/Trinotate.sqlite", 
        annotation_table="24.Generate_annotation_table/{}.trino_anno_rep.tsv".format(config['project_name'])
#    log: "logs/Generate_annotation_table/Generate_annotation_table.log"
    threads: 10
    params:
        conda_activate=config['conda']['non_model_RNA_Seq']['env'],
        PERL5LIB=config['conda']['non_model_RNA_Seq']['perl5lib'],
        program=config['programs_path']['trinotate']['trinotate'],
        threads=30,
        out_dir=lambda w, output: path.dirname(output.sqlitedb),
        database=config['databases']['trinotate']['sqlitedb'],
        o="logs/Generate_annotation_table/Generate_annotation_table.log.o",
        e="logs/Generate_annotation_table/Generate_annotation_table.log.e"
    shell:
        """
        set +u
        {params.conda_activate}
        {params.PERL5LIB}
        set -u

        # Copy Trinotate SQLite Databse to the annotation directory

        cp {params.database} {params.out_dir}
        
        # Intialize database

        {params.program} \
         {output.sqlitedb} \
         init \
         --gene_trans_map {input.gene_trans_map} \
         --transcript_fasta  {input.transcriptome} \
         --transdecoder_pep {input.coding_sequences} \

        # ----- Load reports --------
        # Swissprot blastp

        {params.program} \
         {output.sqlitedb} \
         LOAD_swissprot_blastp  \
         {input.swissprot_blastp}  

        # Swissprot blastx

        {params.program} \
         {output.sqlitedb} \
         LOAD_swissprot_blastx  \
         {input.swissprot_blastx}

        # pfam

        {params.program} \
         {output.sqlitedb} \
         LOAD_pfam  \
         {input.pfam}

        # Ribosomal RNA with enammer

        {params.program} \
         {output.sqlitedb} \
         LOAD_rnammer  \
         {input.rRNA}

        # Create Excel report

        {params.program} \
         {output.sqlitedb} \
         report \
          > {output.annotation_table}
         """


   
