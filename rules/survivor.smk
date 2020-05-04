SAMPLES = config["samples"].keys()

rule survivor:
    input:
        [f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/{sample}/{sample}.vcf"
            for sample in SAMPLES]
    output:
        vcf = f"{OUTDIR}/{{aligner}}/{{caller}}_combined/{{stage}}.vcf",
        fofn = f"{OUTDIR}/{{aligner}}/{{caller}}_combined/{{stage}}.fofn"
    params:
        distance = config["survivor"]["distance"],
        caller_support = config["survivor"]["caller_support"],
        same_type = config["survivor"]["same_type"],
        same_strand = config["survivor"]["same_strand"],
        estimate_distance = config["survivor"]["estimate_distance"],
        minimum_size = config["survivor"]["minimum_size"],
    log:
        out = f"{LOGDIR}/{{aligner}}/{{caller}}/survivor_{{stage}}.out",
        err = f"{LOGDIR}/{{aligner}}/{{caller}}/survivor_{{stage}}.err"
    shell:
        """
        ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
            {params.same_type} {params.same_strand} {params.estimate_distance}  \
            {params.minimum_size} {output.vcf} 1> {log.out} 2> {log.err}
        """

rule survivor_intra_sample:
    input:
        f"{OUTDIR}/{{aligner}}/all_concat/{{sample}}/{{sample}}.vcf"
    output:
        vcf = f"{OUTDIR}/{{aligner}}/all_merged/{{sample}}/{{sample}}.vcf",
        fofn = f"{OUTDIR}/{{aligner}}/all_merged/{{sample}}/{{sample}}.fofn"
    params:
        distance = config["survivor"]["distance"],
        caller_support = 1,
        same_type = -1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        f"{LOGDIR}/{{aligner}}/all_merged/survivor.log"
    shell:
        """
        cat <(ls {input}) > {output.fofn}; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}
        """

# rule survivor_caller:
#     input:
#         [f"{OUTDIR}/{{aligner}}/{{caller}}_reformatted/{sample}/{sample}.vcf"
#             for sample in SAMPLES]
#     output:
#         vcf = f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf",
#         fofn = f"{OUTDIR}/{{aligner}}/{{caller}}_combined/samples.fofn"
#     params:
#         distance = config["survivor"]["distance"],
#         caller_support = 1,
#         same_type = 1,
#         same_strand = -1,
#         estimate_distance = -1,
#         minimum_size = -1,
#     log:
#         f"{LOGDIR}/{{aligner}}/{{caller}}_combined/{{caller}}.log"
#     shell:
#         """
#         cat <(ls {input}) > {output.fofn} ; \
#         SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
#         {params.same_type} {params.same_strand} {params.estimate_distance}  \
#         {params.minimum_size} {output.vcf} 2> {log}
#         """

rule survivor_inter_sample:
    input:
        [f"{OUTDIR}/{{aligner}}/all_merged/{sample}/{sample}.vcf"
               for sample in SAMPLES]
    output:
        vcf = f"{OUTDIR}/{{aligner}}/all_combined/genotypes.vcf",
        fofn = f"{OUTDIR}/{{aligner}}/all_combined/samples.fofn"
    params:
        distance = config["survivor"]["distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
    log:
        f"{LOGDIR}/{{aligner}}/all/survivor.log"
    shell:
        """
        cat <(ls {input}) > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}
        """