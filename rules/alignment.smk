def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule minimap2_align:
    input:
        fq = f"{OUTDIR}/fastq/{{sample}}.fastq.gz",
        genome = config["genome"]
    output:
        f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam"
    threads:
        config["threads"]["per_sample"]
    log:
        f"{LOGDIR}/alignment_sorted/{{sample}}.log"
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
         -R "@RG\\tID:{{sample}}\\tSM:{{sample}}" \
         {input.genome} {input.fq} | \
         samtools sort -@ {threads} -O BAM -o {output} - 2> {log}
        """

rule samtools_sort:
    input:
        f"{OUTDIR}/alignment_unsorted/{{sample}}.bam"
    output:
        f"{OUTDIR}/alignment_sorted/{{sample}}.bam"
    threads:
        config["threads"]["per_sample"]
    log:
        f"{LOGDIR}/samtools_sort/{{sample}}.log"
    shell:
        """
        samtools sort -@ {threads} {input} -O BAM -o {output} 2> {log}
        """

rule samtools_index:
    input:
        f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam"
    output:
        f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam.bai"
    log:
        f"{LOGDIR}/samtools_index/{{sample}}.log"
    shell:
        """
        samtools index {input} 2> {log}
        """

