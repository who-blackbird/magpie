def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

SAMPLES = config["samples"].keys()

rule nanopolish:
    input:
        fq = f"{OUTDIR}/fastq/{{sample}}.fastq.gz",
        bam = f"{OUTDIR}/alignment_sorted/{{sample}}.bam",
        bai = f"{OUTDIR}/alignment_sorted/{{sample}}.bam.bai",
        genome = config["genome"],
        readdb = f"{OUTDIR}/fastq/{{sample}}.fastq.gz.index.readdb" 
    output:
        f"{OUTDIR}/methylation_parallel/{{sample}}/{{range}}"
    threads:
        config["threads"]["per_sample"]
    params:
        os.path.join(workflow.basedir, "scripts/nanopolish_makerange.py")
    log:
        err = f"{LOGDIR}/methylation/{{sample}}.err",
        out = f"{LOGDIR}/methylation/{{sample}}.out" 
    shell:
        """
        nanopolish call-methylation \
            -w {{input.range}} \
            --methylation=cpg \
            --progress \
            -t {threads} \
            -r {input.fq} \
            -b {input.bam} \
            -g {input.genome} > {output} 1> {log.out} 2> {log.err}
        """

rule merge_methylation:
    input:
        range = [f"{OUTDIR}/methylation_parallel/{{sample}}/{range}" for range in RANGES],
        header = path/to/header
    output:
        f"{OUTDIR}/methylation/{{sample}}.tab"
    shell:
        """
        cat <(cat {header}) <(tail -n+2 {input}) > {output}
        """

rule frequency:
    input:
        f"{OUTDIR}/methylation/{{sample}}.tab"
    output:
        f"{OUTDIR}/methylation_frequencies/{{sample}}_frequency.tab"
    params:
        os.path.join(workflow.basedir, "scripts/calculate_methylation_frequency.py")
    log:
        err = f"{LOGDIR}/methylation_freq/{{sample}}.err"
    shell:
        """
        {params} {input} > {output}  2> {log.err}
        """