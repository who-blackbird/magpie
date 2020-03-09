rule alignment_qc:
    input:
        bam = "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    output:
        cov = "{OUTDIR}/{{aligner}}/coverage/{{sample}}.cov",
        stats = "{OUTDIR}/{{aligner}}/coverage/{{sample}}.stats"
    params:
        os.path.join(workflow.basedir, "scripts/coverageStats.py")
    log:
        "{LOGDIR}/{{aligner}}/coverage/{{sample}}.log"
    shell:
        "samtools depth {input.bam} > {output.cov} 2> {log}; \
         python {params} {output.cov} > {output.stats}"
