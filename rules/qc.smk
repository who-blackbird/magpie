rule alignment_qc:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam.bai"
    output:
        cov = temp(f"{OUTDIR}/{{aligner}}/qc/coverage/{{sample}}.cov"),
        stats = f"{OUTDIR}/{{aligner}}/qc/coverage/{{sample}}.stats"
    params:
        os.path.join(workflow.basedir, "scripts/coverage_stats.py")
    log:
        f"{LOGDIR}/{{aligner}}/coverage/{{sample}}.log"
    shell:
        "samtools depth {input.bam} > {output.cov} 2> {log}; \
         python {params} {output.cov} > {output.stats}"

rule nanostat:
    input:
        f"{OUTDIR}/fastq/{{sample}}.fastq.gz"
    output:
        dir = f"{OUTDIR}/{{aligner}}/qc/stats/{{sample}}",
        file = f"{OUTDIR}/{{aligner}}/qc/stats/{{sample}}/{{sample}}.stats"
    log:
        f"{LOGDIR}/{{aligner}}/stats/{{sample}}.log"
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
        f"{OUTDIR}/{{aligner}}/qc/read_length/{{sample}}.txt"
    params:
        os.path.join(workflow.basedir, "scripts/read_length.sh")
    log:
        f"{LOGDIR}/{{aligner}}/read_length/{{sample}}.log"
    shell:
        """
        sh {params} {input} > {output} 2> {log}
        """
