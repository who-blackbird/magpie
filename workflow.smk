import os
import gzip
import sys

#configfile: "config.yaml"

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

include: "rules/fqmerge.smk"
include: "rules/alignment.smk"
include: "rules/svcalling.smk"
include: "rules/qc.smk"
include: "rules/snps.smk"

# Functions #

def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

#CHROMOSOMES = getChr()
CHROMOSOMES = ['1', '4']

# Target rules #

rule structural_variants:
    input:
        expand(f"{OUTDIR}/minimap2/coverage/{{sample}}.stats",
               sample=config["samples"]),
        f"{OUTDIR}/minimap2/sniffles_annotated/genotypes_annot.vcf.gz"


rule snps:
    input:
        expand(f"{OUTDIR}/minimap2/longshot_annotated/all-{{chromosome}}.snps.annot.vcf.gz",
            chromosome=CHROMOSOMES)
