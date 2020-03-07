def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule minimap2_align:
    input:
        fq = f"{OUTDIR}/fastq/{{sample}}.fastq",
        genome = config["genome"]
    output:
        f"{OUTDIR}/minimap2/alignment/{{sample}}.bam"
    threads:
        config["threads"]
    log:
        f"{LOGDIR}/minimap2/{{sample}}.log"
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
         -R "@RG\\tID:{{sample}}\\tSM:{{sample}}" \
         {input.genome} {input.fq} | \
         samtools sort -@ {threads} -o {output} - 2> {log}
        """

#rule ngmlr_align:
    #input:
        #fq = get_samples,
        #genome = config["genome"]
    #output:
        #protected(f"{OUTDIR}/ngmlr/alignment/{{sample}}.bam")
    #threads:
        #config["threads"]
    #log:
        #f"{LOGDIR}/ngmlr/{{sample}}.log"
    #shell:
        #"""
        #zcat {input.fq} | \
        #ngmlr --presets ont -t {threads} -r {input.genome} | \
        #samtools sort -@ {threads} -o {output} - 2> {log}
        #"""

rule samtools_index:
    input:
        f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam"
    output:
        f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    log:
        f"{LOGDIR}/{{aligner}}/samtools_index/{{sample}}.log"
    shell:
        """
        samtools index {input} 2> {log}
        """

