rule vcfanno_svs:
    input:
        f"{OUTDIR}/{{aligner}}/all_merged/genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/all_annotated/svs_vcfanno.vcf"
    log:
        f"{LOGDIR}/{{aligner}}/annotate_all/annotate_vcfanno.err"
    params:
         conf=config["vcfanno_conf_svs"]
    threads:
         config["threads"]["all"]
    shell:
         """
         vcfanno -ends -permissive-overlap -p {threads} {params.conf} {input} > {output} 2> {log}
         """

rule vcf2tab:
    input:
         f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/svs_{{annot}}.vcf"
    output:
          f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/svs_{{annot}}.tab"
    log:
       f"{LOGDIR}/{{aligner}}/vcf2tab/vcf2tab_{{caller}}.err"
    params:
          script=os.path.join(workflow.basedir, "scripts/vcf2tab.py")
    shell:
         """
         python {params.script} {input} > {output}
         """

rule get_overlaps:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/svs_{{annot}}.tab"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/svs_{{annot}}.ovl.tab"
    log:
        f"{LOGDIR}/{{aligner}}/overlaps/overlaps_{{caller}}.err"
    params:
          script=os.path.join(workflow.basedir, "scripts/overlap/main.py"),
          config_file=os.path.join(workflow.basedir, "config/overlap_config.ini"),
          columns_file=os.path.join(workflow.basedir, "config/colnames.yaml")
    shell:
         """
         python {params.script} \
            --config {params.config_file} \
            --columns {params.columns_file} \
            --input {input} \
            --output {output}
         """

         # rule annotsv_svs:
         #     input:
         #         f"{OUTDIR}/{{aligner}}/{{caller}}_combined/missing2ref.vcf"
         #     output:
         #         f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/svs_annotsv.vcf"
         #     log:
         #         err = f"{LOGDIR}/{{aligner}}/annotate_{{caller}}/annotate_annotsv.err"
         #     params:
         #         annotsv = config["annotsv"],
         #         annotsv_path = config["annotsv_path"]
         #     shell:
         #         """
         #         export ANNOTSV={params.annotsv_path}; \
         #         {params.annotsv} -SVinputFile {input} \
         #         -genomeBuild GRCh38 \
         #         -typeOfAnnotation full \
         #         -overlap 70 \
         #         -reciprocal yes \
         #         -outputFile {output} 2> {log.err}
         #         """
