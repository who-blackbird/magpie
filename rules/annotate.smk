rule vcfanno_svs:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_vcfanno.vcf"
    log:
        err = f"{LOGDIR}/{{aligner}}/annotate_svs/annotate_{{caller}}.err"
    params:
        conf = config["vcfanno_conf_svs"]
    threads:
        config["threads"]["all"]
    shell:
        """
        vcfanno -ends -permissive-overlap -p {threads} {params.conf} {input} > {output} 2> {log.err}
        """


rule annotSV:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_annotsv.vcf"
    log:
        err = f"{LOGDIR}/{{aligner}}/annotate_svs/annotate_{{caller}}.err"
    params:
        annotsv = config["annotsv"],
        annotsv_path = config["annotsv_path"]
    shell:
        """
        export ANNOTSV={params.annotsv_path}; \
        {params.annotsv} -SVinputFile {input} \
        -genomeBuild GRCh38 \
        -typeOfAnnotation full \
        -overlap 70 \
        -reciprocal yes \
        -outputFile {output} 2> {log.err}
        """

rule vcf2tab:
    input:
         f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_{{annot}}.vcf"
    output:
          f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_{{annot}}.tab"
    log:
       err = f"{LOGDIR}/{{aligner}}/vcf2tab/vcf2tab_{{caller}}.err"
    params:
        os.path.join(workflow.basedir, "scripts/vcf2tab.py")
    shell:
         """
         python {script} {input} > {output}
         """