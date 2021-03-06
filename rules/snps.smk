SAMPLES = config["samples"].keys()

rule samtools_split:
    input:
        bam = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam",
        bai = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam.bai",
    output:
        bam = f"{OUTDIR}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam",
        bai = f"{OUTDIR}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam.bai"
    log:
        f"{LOGDIR}/samtools_split/{{sample}}/{{sample}}-{{chromosome}}.log"
    shell:
        """
        samtools view -h {input.bam} {wildcards.chromosome} -o {output.bam} 2> {log}
        samtools index {output.bam}
        """

rule longshot_call:
    input:
        bam = f"{OUTDIR}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam",
        bai = f"{OUTDIR}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam.bai",
    output:
        vcf = temp(f"{OUTDIR}/longshot_split/{{sample}}/{{sample}}-{{chromosome}}.snps.vcf"),
    log:
        f"{LOGDIR}/longshot_split/{{sample}}/{{sample}}-{{chromosome}}.snps.log"
    params:
        genome = config["genome"]
    shell:
        """
        longshot --bam {input.bam} --ref {params.genome} --out {output.vcf} 2> {log}
        """

rule bgzip_and_tabix:
    input:
        f"{OUTDIR}/{{caller}}_split/{{sample}}/{{sample}}-{{chromosome}}.snps.vcf"
    output:
        vcf = f"{OUTDIR}/{{caller}}_split/{{sample}}/{{sample}}-{{chromosome}}.snps.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_split/{{sample}}/{{sample}}-{{chromosome}}.snps.vcf.gz.tbi"
    log:
        f"{LOGDIR}/{{caller}}_index/{{sample}}/{{sample}}-{{chromosome}}.index.log"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule merge_vcfs:
    input:
        vcf = [f"{OUTDIR}/{{caller}}_split/{sample}/{sample}-{{chromosome}}.snps.vcf.gz" for sample in SAMPLES],
        idx = [f"{OUTDIR}/{{caller}}_split/{sample}/{sample}-{{chromosome}}.snps.vcf.gz.tbi" for sample in SAMPLES]
    output:
        vcf = f"{OUTDIR}/{{caller}}_merged/all-{{chromosome}}.snps.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_merged/all-{{chromosome}}.snps.vcf.gz.tbi"
    log:
       f"{LOGDIR}/{{caller}}_merge/all-{{chromosome}}.log"
    shell:
        """
        bcftools merge --force-samples {input.vcf} -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule vcfanno:
    input:
        vcf = f"{OUTDIR}/{{caller}}_merged/all-{{chromosome}}.snps.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_merged/all-{{chromosome}}.snps.vcf.gz.tbi"
    output:
        vcf = f"{OUTDIR}/{{caller}}_vcfanno_annotated/all-{{chromosome}}.snps.annot.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_vcfanno_annotated/all-{{chromosome}}.snps.annot.vcf.gz.tbi"
    log:
        f"{LOGDIR}/{{caller}}_{{annot}}annotation/all-{{chromosome}}.annot.log"
    params:
        conf = config["vcfanno_conf_snvs"]
    threads:
        config["threads"]["per_chr"]
    shell:
        """
        vcfanno -p {threads} {params.conf} {input} | \
        bgzip -c > {output} 2> {log.err}
        """

rule vep:
    input:
        vcf = f"{OUTDIR}/{{caller}}_merged/all-{{chromosome}}.snps.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_merged/all-{{chromosome}}.snps.vcf.gz.tbi"
    output:
        vcf = f"{OUTDIR}/{{caller}}_vep_annotated/all-{{chromosome}}.snps.annot.vcf.gz",
        idx = f"{OUTDIR}/{{caller}}_vep_annotated/all-{{chromosome}}.snps.annot.vcf.gz.tbi"
    params:
        os.path.join(workflow.basedir, "scripts/vep_annotation.sh"),
    log:
        f"{LOGDIR}/{{caller}}_vep_annotation/all-{{chromosome}}.annot.log"
    shell:
        """
        {params} {input.vcf} {output.vcf} {wildcards.chromosome}
        """
