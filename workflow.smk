import os
import gzip
import sys

#configfile: "config.yaml"

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

include: "rules/fqmerge.smk"
include: "rules/alignment.smk"
include: "rules/svcalling.smk"
include: "rules/survivor.smk"
include: "rules/annotate.smk"
include: "rules/qc.smk"
include: "rules/snps.smk"
include: "rules/vcf.smk"

# Functions #

def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

CHROMOSOMES = getChr()

# Target rules #

rule fast:
    input:
        expand(f"{OUTDIR}/minimap2/qc/coverage/{{sample}}.stats",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/qc/stats/{{sample}}/{{sample}}.stats",
               sample=config["samples"]),
        # f"{OUTDIR}/minimap2/sniffles_annotated/svs_vcfanno.ovl.tab"

rule precise:
    input:
        expand(f"{OUTDIR}/minimap2/qc/coverage/{{sample}}.stats",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/qc/stats/{{sample}}/{{sample}}.stats",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/qc/read_length/{{sample}}.txt",
               sample=config["samples"]),
        expand(f"{OUTDIR}/minimap2/all_filtered/{{sample}}/{{sample}}.vcf",
               sample=config["samples"]),
        # f"{OUTDIR}/minimap2/all_annotated/svs_vcfanno.ovl.tab"

rule snps:
    input:
        expand(f"{OUTDIR}/minimap2/longshot_split/{{sample}}/{{sample}}-{{chromosome}}.snps.vcf.gz",
            sample=config["samples"], chromosome=CHROMOSOMES)
        #expand(f"{OUTDIR}/minimap2/longshot_split/{{sample}}/{{sample}}-{{chromosome}}.snps.vcf",
            # sample=config["samples"], chromosome=CHROMOSOMES)
        #expand(f"{OUTDIR}/minimap2/longshot_vep_annotated/all-{{chromosome}}.snps.annot.vcf.gz",
            #chromosome=CHROMOSOMES)
