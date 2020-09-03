rule vcfanno_svs:
    input:
        f"{OUTDIR}/all_merged/genotypes.vcf"
    output:
        f"{OUTDIR}/all_annotated/svs_vcfanno.vcf"
    log:
        f"{LOGDIR}/annotate_all/annotate_vcfanno.err"
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
         f"{OUTDIR}/{{caller}}_annotated/svs_{{annot}}.vcf"
    output:
          f"{OUTDIR}/{{caller}}_annotated/svs_{{annot}}.tab"
    log:
       f"{LOGDIR}/vcf2tab/vcf2tab_{{caller}}.err"
    params:
          script=os.path.join(workflow.basedir, "scripts/vcf2tab.py")
    shell:
         """
         python {params.script} {input} > {output}
         """

rule get_overlaps:
    input:
        f"{OUTDIR}/{{caller}}_annotated/svs_{{annot}}.tab"
    output:
        f"{OUTDIR}/{{caller}}_annotated/svs_{{annot}}.ovl.tab"
    log:
        err=f"{LOGDIR}/overlaps/overlaps_{{caller}}.err",
        out=f"{LOGDIR}/overlaps/overlaps_{{caller}}.out"
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
            --output {output} \
            --verbose 1> {log.out} 2> {log.err}
         """
