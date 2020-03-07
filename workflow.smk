import os
import gzip
import sys

configfile: "config.yaml"

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

include: "rules/fqmerge.smk"
include: "rules/alignment.smk"
include: "rules/svcalling.smk"
include: "rules/annotate.smk"
include: "rules/qc.smk"
include: "rules/snps.smk"

# Target rules #

rule fast:
    input:
        expand(f"{OUTDIR}/minimap2/coverage/{{sample}}.stats",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/sniffles_SVs/{{sample}}.vcf",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/sniffles_annotate/{{sample}}_annot.vcf.tsv",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/longshot/{{sample}}.merged.vcf.gz",
               sample=config["samples"])
