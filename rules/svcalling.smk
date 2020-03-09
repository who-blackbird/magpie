rule sniffles_call:
    input:
        bam = "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    output:
        "{OUTDIR}/{{aligner}}/sniffles_svs/{{sample}}.vcf"
    params:
        se = config["sniffles_se"],
    threads: 
        config["threads"]
    log:
        "{LOGDIR}/{{aligner}}/sniffles_svs/{{sample}}.log"
    shell:
        "sniffles -s {params.se} --mapped_reads {input.bam} --vcf {output} --threads {threads} 2> {log}"
