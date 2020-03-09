def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule fq_merge:
    input:
        get_samples
    output:
        "{OUTDIR}/fastq/{{sample}}.fastq"
    shell:
        """
        cat {input} | xargs cat > {output}
        """
