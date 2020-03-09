def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

CHROMOSOMES = getChr()

rule samtools_split:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai",
    output:
        bam = temp(f"{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam"),
        bai = temp(f"{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam.bai")
    log:
        f"{LOGDIR}/{{aligner}}/samtools_split/{{sample}}-{{chromosome}}.log"
    shell:
        """
        samtools view -h {input.bam} {wildcards.chromosome} -o {output.bam} 2> {log}
        samtools index {output.bam}
        """

rule longshot_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_split/{{sample}}-{{chromosome}}.bam.bai",
    output:
        vcf   = temp(f"{OUTDIR}/{{aligner}}/longshot_split/{{sample}}-{{chromosome}}.snps.vcf"),
    log:
        f"{LOGDIR}/{{aligner}}/longshot_split/{{sample}}-{{chromosome}}.snps.log"
    params:
        genome = config["genome"]
    shell:
        """
        longshot --bam {input.bam} --ref {params.genome} --out {output.vcf} 2> {log}
        """

rule bgzip_and_tabix:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{{chromosome}}.snps.vcf"
    output:
        vcf = temp(f"{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{{chromosome}}.snps.vcf.gz"),
        idx = temp(f"{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{{chromosome}}.snps.vcf.gz.tbi")
    log:
        f"{LOGDIR}/{{aligner}}/{{caller}}_index/{{sample}}-{{chromosome}}.index.log"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule cat_vcfs:
    input:
        vcf = [f"{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{chromosome}.snps.vcf.gz" for chromosome in CHROMOSOMES],
        idx = [f"{OUTDIR}/{{aligner}}/{{caller}}_split/{{sample}}-{chromosome}.snps.vcf.gz.tbi" for chromosome in CHROMOSOMES]
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}/{{sample}}.merged.vcf.gz"
    log:
        f"{LOGDIR}/{{aligner}}/{{caller}}/bcftools-concat/{{sample}}.merged.log"
    shell:
        """
        bcftools concat {input.vcf} -O z -o {output}
        """
