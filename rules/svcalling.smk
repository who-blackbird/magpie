SAMPLES = config["samples"].keys()

rule sniffles_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment_sorted/{{sample}}.bam.bai"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_calls/{{sample}}.vcf"
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
        [f"{OUTDIR}/{{aligner}}/{{caller}}_{{stage}}/{sample}.vcf" for sample in SAMPLES]
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
        f"{OUTDIR}/{{aligner}}/sniffles_genotypes/{{sample}}.vcf"
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
        f"{OUTDIR}/{{aligner}}/nanosv_genotypes/{{sample}}/{{sample}}-{{chromosome}}.vcf"
    threads:
        2
    params:
        bed = config["annotbed"]
    log:
        out = f"{LOGDIR}/{{aligner}}/nanosv_genotypes/{{sample}}.out",
        err = f"{LOGDIR}/{{aligner}}/nanosv_genotypes/{{sample}}.err"
    shell:
        """
       NanoSV -s sambamba \
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
        # f"{OUTDIR}/{{aligner}}/svim_genotypes/{{sample}}/{{sample,[A-Za-z0-9]+}}.vcf"
        # f"{OUTDIR}/{{aligner}}/svim_genotypes/{{sample}}/{{sample}}.vcf"
    log:
        f"{LOGDIR}/{{aligner}}/svim_genotype/{{sample}}.filter.log"
    shell:
        """
        cat {input} | \
        awk '{{ if($1 ~ /^#/) {{ print $0 }} \
        else {{ if($6>10) {{ print $0 }} }} }}' > {output}
        """

rule missing2ref:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_genotypes/{{sample}}.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_missing2ref/{{sample}}.vcf"
    log:
        f"{LOGDIR}/{{aligner}}/{{caller}}_missing2ref/{{sample}}.err"
    shell:
        """
        bcftools +missing2ref {input} > {output} 2> {log}
        """