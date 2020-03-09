def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule minimap2_align:
    input:
        fq = "{OUTDIR}/fastq/{{sample}}.fastq",
        genome = config["genome"]
    output:
        "{OUTDIR}/minimap2/alignment/{{sample}}.bam"
    threads:
        config["threads"]
    log:
        "{LOGDIR}/minimap2/{{sample}}.log"
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
        #protected("{OUTDIR}/ngmlr/alignment/{{sample}}.bam")
    #threads:
        #config["threads"]
    #log:
        #"{LOGDIR}/ngmlr/{{sample}}.log"
    #shell:
        #"""
        #zcat {input.fq} | \
        #ngmlr --presets ont -t {threads} -r {input.genome} | \
        #samtools sort -@ {threads} -o {output} - 2> {log}
        #"""

rule samtools_index:
    input:
        "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam"
    output:
        "{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    log:
        "{LOGDIR}/{{aligner}}/samtools_index/{{sample}}.log"
    shell:
        """
        samtools index {input} 2> {log}
        """

