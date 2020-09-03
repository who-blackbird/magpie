def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

CHROMOSOMES = getChr()

rule sniffles_call:
    input:
        bam = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam",
        bai = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam.bai"
    output:
        f"{OUTDIR}/sniffles_genotypes/{{sample}}/{{sample}}.vcf"
    params:
        se = config["sniffles_se"],
    threads:
        # config["threads"]["per_sample"]
        1
    log:
        out = f"{LOGDIR}/sniffles_genotypes/{{sample}}.out",
        err = f"{LOGDIR}/sniffles_genotypes/{{sample}}.err"
    shell:
        """
        sniffles -s {params.se} --mapped_reads {input.bam} --vcf {output} --threads {threads} 1> {log.out} 2> {log.err}
        """

rule nanosv_call:
    input:
        bam = f"{OUTDIR}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam",
        bai = f"{OUTDIR}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam.bai"
    output:
        f"{OUTDIR}/nanosv_genotypes_split/{{sample}}/{{sample}}-{{chromosome}}.vcf"
    threads:
        1
    params:
        bed = config["annotbed"],
        sambamba = config["sambamba_path"]
    log:
        f"{LOGDIR}/nanosv_genotypes_split/{{sample}}-{{chromosome}}.log"
    shell:
        """
        NanoSV -s {params.sambamba} \
        -t {threads} \
        -b {params.bed} \
        -o {output} {input.bam} 2> {log} 
        """

rule svim_call:
    input:
        bam = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam",
        bai = f"{OUTDIR}/alignment_sorted/{{sample}}/{{sample}}.bam.bai",
        genome = config["genome"]
    output:
        f"{OUTDIR}/svim_calls/{{sample}}/final_results.vcf"
    params:
        outdir=f"{OUTDIR}/svim_calls/{{sample}}"
    log:
        f"{LOGDIR}/svim_calls/{{sample}}.log"
    shell:
        """
        svim alignment --sample {wildcards.sample} \
        {params.outdir}/ {input.bam} {input.genome} 2> {log}
        """

rule filter_svim:
    input:
        f"{OUTDIR}/svim_calls/{{sample}}/final_results.vcf"
    output:
        f"{OUTDIR}/svim_genotypes/{{sample}}/{{sample}}.vcf"
    log:
        f"{LOGDIR}/svim_genotype/{{sample}}.filter.log"
    shell:
        """
        cat {input} | \
        awk '{{ if($1 ~ /^#/) {{ print $0 }} \
        else {{ if($6>10) {{ print $0 }} }} }}' > {output}
        """
