def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

CHROMOSOMES = getChr()

rule samtools_split:
    input:
        bam = "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai",
    output:
        bam = temp("{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam"),
        bai = temp("{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam.bai")
    log:
        "{LOGDIR}/{{aligner}}/samtools_split/{{sample}}-{{chromosome}}.log"
    shell:
        """
        samtools view -h {input.bam} {wildcards.chromosome} -o {output.bam} 2> {log}
        samtools index {output.bam}
        """

rule longshot_call:
    input:
        bam = "{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam",
        bai = "{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam.bai",
    output:
        vcf   = temp("{OUTDIR}/{{aligner}}/longshot_split/{{sample}}-{{chromosome}}.snps.vcf"),
    log:
        "{LOGDIR}/{{aligner}}/longshot_split/{{sample}}-{{chromosome}}.snps.log"
    params:
        genome = config["genome"]
    shell:
        """
        longshot --bam {input.bam} --ref {params.genome} --out {output.vcf} 2> {log}
        """

rule bgzip_and_tabix:
    input:
        "{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{{chromosome}}.snps.vcf"
    output:
        vcf = temp("{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{{chromosome}}.snps.vcf.gz"),
        idx = temp("{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{{chromosome}}.snps.vcf.gz.tbi")
    log:
        "{LOGDIR}/{{aligner}}/{{caller}}_index/{{sample}}-{{chromosome}}.index.log"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule cat_vcfs:
    input:
        vcf = ["{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{chromosome}.snps.vcf.gz" for chromosome in CHROMOSOMES],
        idx = ["{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{chromosome}.snps.vcf.gz.tbi" for chromosome in CHROMOSOMES]
    output:
        "{OUTDIR}/{{aligner}}/{{caller}}/{{sample}}.merged.vcf.gz"
    log:
        "{LOGDIR}/{{aligner}}/{{caller}}/bcftools-concat/{{sample}}.merged.log"
    shell:
        """
        bcftools concat {input.vcf} -O z -o {output}
        """
