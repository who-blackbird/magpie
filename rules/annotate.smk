rule annotate_vcf:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/{{sample}}.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/{{sample}}_annot.vcf.tsv"
    log:
        f"{LOGDIR}/{{aligner}}/annotate_vcf/annotate_{{caller}}.log"
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
        -outputFile {output}
        """
