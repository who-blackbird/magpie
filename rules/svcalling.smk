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
        f"{LOGDIR}/{{aligner}}/sniffles_calls/{{sample}}.log"
    shell:
        """
        sniffles -s {params.se} --mapped_reads {input.bam} --vcf {output} --threads {threads} 2> {log}
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
        f"{LOGDIR}/{{aligner}}/{{caller}}/survivor_{{stage}}.log"
    shell:
        """
        ls {input} > {output.fofn} ; \
        SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
            {params.same_type} {params.same_strand} {params.estimate_distance}  \
            {params.minimum_size} {output.vcf} 2> {log}
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
        f"{LOGDIR}/{{aligner}}/sniffles_genotypes/{{sample}}.log"
    shell:
        """
        sniffles --mapped_reads {input.bam} \
            --vcf {output} \
            --threads {threads} \
            --report_seq \
            --cluster \
            --Ivcf {input.ivcf} 2> {log}
        """

rule annotate_vcf:
    input:
        f"{OUTDIR}/{{aligner}}/{{caller}}_combined/genotypes.vcf"
    output:
        f"{OUTDIR}/{{aligner}}/{{caller}}_annotated/genotypes_annot.tsv"
    log:
        f"{LOGDIR}/{{aligner}}/annotate_vcf/annotate_{{caller}}.log"
    params:
        annotsv = config["annotsv"],
        annotsv_path = config["annotsv_path"]
    shell:
        """
        export ANNOTSV={params.annotsv_path}; \
        {params.annotsv} -SVinputFile {input} \
            -genomeBuild GRCh38 \
            -typeOfAnnotation full \
            -overlap 70 \
            -reciprocal yes \
            -outputFile {output}
        """
