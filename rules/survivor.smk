SAMPLES = config["samples"].keys()

rule survivor_intra_sample:
    input:
        f"{OUTDIR}/intrasample_concat/{{sample}}/{{sample}}.vcf"
    output:
        vcf = f"{OUTDIR}/intrasample_merged/{{sample}}/{{sample}}.vcf",
        fofn = f"{OUTDIR}/intrasample_merged/{{sample}}/{{sample}}.fofn"
    params:
        SURVIVOR = config["survivor"]["src"],
        distance = config["survivor"]["distance"],
        caller_support = 1,
        same_type = -1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        f"{LOGDIR}/intrasample_merged/survivor.log"
    shell:
        """
        cat <(ls {input}) > {output.fofn}; \
        {params.SURVIVOR} merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}
        """

rule survivor_inter_sample:
    input:
        [f"{OUTDIR}/intrasample_changetype/{sample}/{sample}.vcf"
               for sample in SAMPLES]
    output:
        vcf = f"{OUTDIR}/all_merged/genotypes.vcf",
        fofn = f"{OUTDIR}/all_merged/samples.fofn"
    params:
        SURVIVOR = config["survivor"]["src"],
        distance = config["survivor"]["distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        f"{LOGDIR}/all_merge/survivor.log"
    shell:
        """
        cat <(ls {input}) > {output.fofn} ; \
        {params.SURVIVOR} merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}
        """
