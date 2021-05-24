
rule Humann2_normalize_per_sample:
    input:
        gene_families="10.Humann2_classify/{sample}/{sample}_genefamilies.tsv",
        pathabundance="10.Humann2_classify/{sample}/{sample}_pathabundance.tsv"
    output:
        gene_families="10.Humann2_classify/{sample}/{sample}_genefamilies_relab.tsv",
        pathabundance="10.Humann2_classify/{sample}/{sample}_pathabundance_relab.tsv"
    params:
        conda_activate=config['conda']['metagenomics']['env'],
        program=config['programs_path']['humann2']['renorm'],
        o="logs/Humann2_normalize/{sample}/{sample}.log.o",
        e="logs/Humann2_normalize/{sample}/{sample}.log.e"
    threads: 1
#    log: "logs/Humann2_normalize/{sample}/{sample}.log"
    shell:
        """
        set +u
        {params.conda_activate}
        set -u
        {params.program} \
         -i {input.gene_families} \
         -o {output.gene_families} \
         --units relab

        {params.program} \
         -i {input.pathabundance} \
         -o {output.pathabundance} \
         --units relab
        """

