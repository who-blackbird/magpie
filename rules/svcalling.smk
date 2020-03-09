rule sniffles_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_{{stage}}/{{sample}}.vcf"
    params:
        se = config["sniffles_se"],
    threads: 
        config["threads"]
    log:
        f"{LOGDIR}/{{aligner}}/sniffles_{{stage}}/{{sample}}.log"
    shell:
        "sniffles -s {params.se} --mapped_reads {input.bam} --vcf {output} --threads {threads} 2> {log}"
