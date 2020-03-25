SAMPLES = config["samples"].keys()

rule sniffles_call:
    input:
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        bai = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam.bai"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_calls/{{sample}}.vcf"
    params:
        se = config["sniffles_se"],
    threads: 
        config["threads"]
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
        distance = config["survivor_distance"],
        caller_support = 1,
        same_type = 1,
        same_strand = -1,
        estimate_distance = -1,
        minimum_size = -1,
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
        bam = f"{OUTDIR}/{{aligner}}/alignment/{{sample}}.bam",
        ivcf = f"{OUTDIR}/{{aligner}}/sniffles_combined/calls.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/sniffles_genotypes/{{sample}}.vcf"
    threads: 
        config["threads"]
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

rule annotate_vcf:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_annot.vcf.gz"
    log:
        err = f"{LOGDIR}/{{aligner}}/annotate_vcf/annotate_{{caller}}.err"
    params:
        conf = config["vcfanno_conf_svs"]
    threads:
        config["threads"]
    shell:
        """
        vcfanno -ends -permissive-overlap -p {threads} {params.conf} {input} | \
        bgzip -c > {output} 2> {log.err}
        """
