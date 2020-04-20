SAMPLES = config["samples"].keys()

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

rule survivor:
    input:
        [f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/{sample}/{sample}.vcf" for sample in SAMPLES]
    output:
        vcf = f"{OUTDIR}/{{aligner}}/{{caller}}_combined/{{stage}}.vcf",
        fofn = temp(f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/samples.fofn")
    params:
        distance = config["survivor"]["distance"],
        caller_support = config["survivor"]["caller_support"],
        same_type = config["survivor"]["same_type"],
        same_strand = config["survivor"]["same_strand"],
        estimate_distance = config["survivor"]["estimate_distance"],
        minimum_size = config["survivor"]["minimum_size"],
    log:
        out = f"{LOGDIR}/{{aligner}}/{{caller}}/survivor_{{stage}}.out",
        err = f"{LOGDIR}/{{aligner}}/{{caller}}/survivor_{{stage}}.err"
    shell:
        """
        ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
            {params.same_type} {params.same_strand} {params.estimate_distance}  \
            {params.minimum_size} {output.vcf} 1> {log.out} 2> {log.err}
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

rule cat_vcfs:
    input:
        [f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes_split/{{sample}}/{{sample}}-{chromosome}.vcf" for chromosome in CHROMOSOMES]
    output:
        unsorted = temp(f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes/{{sample}}/{{sample}}.unsorted.vcf"),
        sorted = f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes/{{sample}}/{{sample}}.vcf"
    log:
        concat = f"{LOGDIR}/{{aligner}}/{{caller}}_concat/{{sample}}.log",
        sort = f"{LOGDIR}/{{aligner}}/{{caller}}_sort/{{sample}}.log"
    shell:
        """
        bcftools concat {input} -o {output.unsorted} 2> {log.concat}; \
        bcftools sort {output.unsorted} -o {output.sorted} 2 > {log.sort}
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

rule vcf_ref:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes/{{sample}}/{{sample}}.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_reformatted/{{sample}}/{{sample}}.vcf"
    log:
        f"{LOGDIR}/{{aligner}}/{{caller}}_reformatted/{{sample}}.err"
    shell:
        """
        bcftools reheader -s <(echo {wildcards.sample}) {input} | \
        bcftools +missing2ref - > {output} 2> {log}
        """

rule survivor_all:
    input:
        [f"{OUTDIR}/{{aligner}}/{caller}_genotypes/{sample}/{sample}.vcf"
               for sample in config["samples"]
               for caller in ["sniffles", "svim", "nanosv"]]
    output:
        vcf = f"{OUTDIR}/{{aligner}}/all_combined/genotypes.vcf",
        fofn = f"{OUTDIR}/{{aligner}}/all_combined/samples.fofn"
    params:
        svs_short = config["svs_short_fof"],
        distance = config["survivor"]["distance"],
        caller_support = config["survivor"]["caller_support"],
        same_type = config["survivor"]["same_type"],
        same_strand = config["survivor"]["same_strand"],
        estimate_distance = config["survivor"]["estimate_distance"],
        minimum_size = config["survivor"]["minimum_size"],
    log:
        f"{LOGDIR}/{{aligner}}/all/survivor.log"
    shell:
        "cat <(ls {input}) <(cat {params.svs_short}) > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
        {params.same_type} {params.same_strand} {params.estimate_distance}  \
        {params.minimum_size} {output.vcf} 2> {log}"