rule vcfanno_svs:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_vcfanno.vcf.gz"
    log:
        err = f"{LOGDIR}/{{aligner}}/annotate_svs/annotate_{{caller}}.err"
    params:
        conf = config["vcfanno_conf_svs"]
    threads:
        config["threads"]["all"]
    shell:
        """
        vcfanno -ends -permissive-overlap -p {threads} {params.conf} {input} | \
        bgzip -c > {output} 2> {log.err}
        """


rule annotSV:
    input:
        "{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        "{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_annotsv.vcf"
    log:
        err = "{LOGDIR}/{{aligner}}/annotate_svs/annotate_{{caller}}.err"
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