rule cat_vcfs:
    input:
        [f"{OUTDIR}/{{caller}}_genotypes_split/{{sample}}/{{sample}}-{chromosome}.vcf" for chromosome in CHROMOSOMES]
    output:
        unsorted = temp(f"{OUTDIR}/{{caller}}_genotypes/{{sample}}/{{sample}}.unsorted.vcf"),
        sorted = f"{OUTDIR}/{{caller}}_genotypes/{{sample}}/{{sample}}.vcf"
    log:
        concat = f"{LOGDIR}/{{caller}}_concat/{{sample}}.log",
        sort = f"{LOGDIR}/{{caller}}_sort/{{sample}}.log"
    shell:
        """
        bcftools concat {input} -o {output.unsorted} 2> {log.concat}; \
        bcftools sort {output.unsorted} -o {output.sorted} 2 > {log.sort}
        """

rule vcf_reheader:
    input:
        f"{OUTDIR}/{{caller}}_genotypes/{{sample}}/{{sample}}.vcf"
    output:
        f"{OUTDIR}/{{caller}}_reheader/{{sample}}/{{sample}}.vcf"
    log:
        f"{LOGDIR}/{{caller}}_reheader/{{sample}}.err"
    shell:
        """
        bcftools reheader -s <(echo {wildcards.sample}) {input} > {output}
        """

rule vcf_ref:
    input:
        f"{OUTDIR}/{{caller}}_reheader/{{sample}}/{{sample}}.vcf"
    output:
        vcf = f"{OUTDIR}/{{caller}}_reformatted/{{sample}}/{{sample}}.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_reformatted/{{sample}}/{{sample}}.vcf.gz.tbi"
    params:
        os.path.join(workflow.basedir, "scripts/change_sv_id.py")
    log:
        f"{LOGDIR}/{{caller}}_reformatted/{{sample}}.err"
    shell:
        """
        python {params} --vcf {input} --caller {wildcards.caller} | bcftools sort - -O z -o {output.vcf}; \
        tabix -p vcf {output.vcf} 2> {log}
        """


rule vcf_concat:
    input:
        [f"{OUTDIR}/{caller}_reformatted/{{sample}}/{{sample}}.vcf.gz"
         for caller in ["sniffles", "svim", "nanosv"]]
    output:
        f"{OUTDIR}/intrasample_concat/{{sample}}/{{sample}}.vcf"
    params:
        svs_short = config["svs_short_fof"]
    log:
        f"{LOGDIR}/intrasample_concat/{{sample}}.err"
    shell:
        """
        bcftools concat -a {input} `grep {wildcards.sample} {params.svs_short} | xargs` | \
        bcftools sort - -o {output}
        """

# rule vcf_blacklist:
#     input:
#         f"{OUTDIR}/intrasample_merged/{{sample}}/{{sample}}.vcf"
#     output:
#         f"{OUTDIR}/intrasample_filtered/{{sample}}/{{sample}}.vcf"
#     params:
#         blacklist = config["blacklist"]
#     log:
#         f"{LOGDIR}/intrasample_filter/{{sample}}.err"
#     shell:
#         """
#         cat <(grep ^# {input}) <(bedtools intersect -a {input} -b {params} -v -f 0.5) > {output}
#         """

rule vcf_change_type:
    input:
        f"{OUTDIR}/intrasample_merged/{{sample}}/{{sample}}.vcf"
    output:
        f"{OUTDIR}/intrasample_changetype/{{sample}}/{{sample}}.vcf"
    params:
        os.path.join(workflow.basedir, "scripts/change_sv_type.py")
    log:
        f"{LOGDIR}/intrasample_changetype/{{sample}}.err"
    shell:
        """
        python {params} --vcf {input} --out {output} 2> {log}
        """
