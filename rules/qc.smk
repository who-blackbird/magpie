rule alignment_qc:
    input:
        bam = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam",
        bai = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam.bai"
    output:
        cov = temp(f"{OUTDIR}/qc/coverage/{{sample}}.cov"),
        stats = f"{OUTDIR}/qc/coverage/{{sample}}.stats"
    params:
        os.path.join(workflow.basedir, "scripts/coverage_stats.py")
    log:
        f"{LOGDIR}/coverage/{{sample}}.log"
    shell:
        "samtools depth {input.bam} > {output.cov} 2> {log}; \
         python {params} {output.cov} > {output.stats}"

rule nanostat:
    input:
        f"{OUTDIR}/fastq/{{sample}}.fastq.gz"
    output:
        dir = f"{OUTDIR}/qc/stats/{{sample}}",
        file = f"{OUTDIR}/qc/stats/{{sample}}/{{sample}}.stats"
    log:
        f"{LOGDIR}/stats/{{sample}}.log"
    threads:
        config["threads"]["per_sample"]
    shell:
        """
        NanoStat --fastq {input} \
                --outdir {output.dir} \
                --prefix {wildcards.sample} \
                --name {output.file} \
                --threads {threads} 2> {log}
        """

rule read_length:
    input:
        f"{OUTDIR}/fastq/{{sample}}.fastq.gz"
    output:
        f"{OUTDIR}/qc/read_length/{{sample}}.txt"
    params:
        os.path.join(workflow.basedir, "scripts/read_length.sh")
    log:
        f"{LOGDIR}/read_length/{{sample}}.log"
    shell:
        """
        sh {params} {input} > {output} 2> {log}
        """
