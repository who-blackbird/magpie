rule alignment_qc:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam.bai"
    output:
        cov = temp(f"{OUTDIR}/{{aligner}}/coverage/{{sample}}.cov"),
        stats = f"{OUTDIR}/{{aligner}}/coverage/{{sample}}.stats"
    params:
        os.path.join(workflow.basedir, "scripts/coverage_stats.py")
    log:
        "{LOGDIR}/{{aligner}}/coverage/{{sample}}.log"
    shell:
        "samtools depth {input.bam} > {output.cov} 2> {log}; \
         python {params} {output.cov} > {output.stats}"
