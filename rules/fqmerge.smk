def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule fq_merge:
    input:
        get_samples
    output:
        f"{OUTDIR}/fastq/{{sample}}.fastq.gz"
    shell:
        """
        cat {input} | xargs cat | gzip > {output}
        """

rule fq_f5_annotation:
    input:
        fq = f"{OUTDIR}/fastq/{{sample}}.fastq.gz",
        f5 = f"{OUTDIR}/fast5/{{sample}}",
    output:
        f"{OUTDIR}/fastq/{{sample}}.fastq.gz.index.readdb"
    shell:
        """
        nanopolish index -d {input.f5} {input.fq}
        """