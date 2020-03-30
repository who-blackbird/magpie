def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule minimap2_align:
    input:
        fq = f"{OUTDIR}/fastq/{{sample}}.fastq.gz",
        genome = config["genome"]
    output:
        f"{OUTDIR}/minimap2/alignment_sorted/{{sample}}.bam"
    threads:
        config["threads"]["per_sample"]
    log:
        f"{LOGDIR}/minimap2/alignment_sorted/{{sample}}.log"
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
         -R "@RG\\tID:{{sample}}\\tSM:{{sample}}" \
         {input.genome} {input.fq} | \
         samtools sort -@ {threads} -O BAM -o {output} - 2> {log}
        """

rule ngmlr_align:
    input:
        fq = f"{OUTDIR}/fastq/{{sample}}.fastq.gz",
        genome = config["genome"]
    output:
        temp(f"{OUTDIR}/ngmlr/alignment_unsorted/{{sample}}.bam")
    threads:
        config["threads"]["per_sample"]
    log:
        f"{LOGDIR}/ngmlr/alignment_unsorted/{{sample}}.log"
    shell:
        """
        ngmlr -x ont -t {threads} -r {input.genome} -q {input.fq} | \
        samtools view - -O BAM -o {output} 2> {log}
        """

rule samtools_sort:
    input:
        f"{OUTDIR}/{{aligner}}/alignment_unsorted/{{sample}}.bam"
    output:
        f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam"
    threads:
        config["threads"]["per_sample"]
    log:
        f"{LOGDIR}/{{aligner}}/samtools_sort/{{sample}}.log"
    shell:
        """
        samtools sort -@ {threads} {input} -O BAM -o {output} 2> {log}
        """

rule samtools_index:
    input:
        f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam"
    output:
        f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam.bai"
    log:
        f"{LOGDIR}/{{aligner}}/samtools_index/{{sample}}.log"
    shell:
        """
        samtools index {input} 2> {log}
        """

