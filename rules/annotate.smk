rule annotate_vcf:
    input:
        "{OUTDIR}/{{aligner}}/{{caller}}_SVs/{{sample}}.vcf"
    output:
        "{OUTDIR}/{{aligner}}/{{caller}}_annotate/{{sample}}_annot.vcf.tsv"
    log:
        "{LOGDIR}/{{aligner}}/annotate_vcf/annotate_{{caller}}.log"
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
