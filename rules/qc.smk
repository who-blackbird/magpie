rule alignment_qc:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    output:
        cov = f"{OUTDIR}/{{aligner}}/coverage/{{sample}}.cov",
        stats = f"{OUTDIR}/{{aligner}}/coverage/{{sample}}.stats"
    params:
        os.path.join(workflow.basedir, "scripts/coverageStats.py")
    log:
        f"{LOGDIR}/{{aligner}}/coverage/{{sample}}.log"
    shell:
        "samtools depth {input.bam} > {output.cov} 2> {log}; \
         python {params} {output.cov} > {output.stats}"
