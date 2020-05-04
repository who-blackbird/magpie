def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

CHROMOSOMES = getChr()

rule sniffles_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam.bai"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_calls/{{sample}}/{{sample}}.vcf"
    params:
        se = config["sniffles_se"],
    threads:
        config["threads"]["per_sample"]
    log:
        out = f"{LOGDIR}/{{aligner}}/sniffles_calls/{{sample}}.out",
        err = f"{LOGDIR}/{{aligner}}/sniffles_calls/{{sample}}.err"
    shell:
        """
        sniffles -s {params.se} --mapped_reads {input.bam} --vcf {output} --threads {threads} 1> {log.out} 2> {log.err}
        """

rule sniffles_genotype:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam",
        ivcf = f"{OUTDIR}/{{aligner}}/sniffles_combined/calls.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_genotypes/{{sample}}/{{sample}}.vcf"
    threads:
        config["threads"]["per_sample"]
    log:
        out = f"{LOGDIR}/{{aligner}}/sniffles_genotypes/{{sample}}.out",
        err = f"{LOGDIR}/{{aligner}}/sniffles_genotypes/{{sample}}.err"
    shell:
        """
        sniffles --mapped_reads {input.bam} \
            --vcf {output} \
            --threads {threads} \
            --report_seq \
            --cluster \
            --Ivcf {input.ivcf} 1> {log.out} 2> {log.err}
        """

rule nanosv_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_split/{{sample}}/{{sample}}-{{chromosome}}.bam.bai"
    output:
        f"{OUTDIR}/{{aligner}}/nanosv_genotypes_split/{{sample}}/{{sample}}-{{chromosome}}.vcf"
    threads:
        1
    params:
        bed = config["annotbed"],
        sambamba = config["sambamba_path"]
    log:
        f"{LOGDIR}/{{aligner}}/nanosv_genotypes_split/{{sample}}-{{chromosome}}.log"
    shell:
        """
        NanoSV -s {params.sambamba} \
        -t {threads} \
        -b {params.bed} \
        -o {output} {input.bam} 2> {log} 
        """

rule svim_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam.bai",
        genome = config["genome"]
    output:
        f"{OUTDIR}/{{aligner}}/svim_calls/{{sample}}/final_results.vcf"
    params:
        outdir=f"{OUTDIR}/{{aligner}}/svim_calls/{{sample}}"
    log:
        f"{LOGDIR}/{{aligner}}/svim_calls/{{sample}}.log"
    shell:
        """
        svim alignment --sample {wildcards.sample} \
        {params.outdir}/ {input.bam} {input.genome} 2> {log}
        """

rule filter_svim:
    input:
        f"{OUTDIR}/{{aligner}}/svim_calls/{{sample}}/final_results.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/svim_genotypes/{{sample}}/{{sample}}.vcf"
    log:
        f"{LOGDIR}/{{aligner}}/svim_genotype/{{sample}}.filter.log"
    shell:
        """
        cat {input} | \
        awk '{{ if($1 ~ /^#/) {{ print $0 }} \
        else {{ if($6>10) {{ print $0 }} }} }}' > {output}
        """