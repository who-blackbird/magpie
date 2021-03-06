## Import libraries
import os
import gzip
import sys

##Directories
OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

##Include rules
include: "rules/fqmerge.smk"
include: "rules/alignment.smk"
include: "rules/svcalling.smk"
include: "rules/survivor.smk"
include: "rules/annotate.smk"
include: "rules/qc.smk"
include: "rules/snps.smk"
include: "rules/vcf.smk"
# include: "rules/methylation.smk" ##In progr

##Functions

def getChr():
    return list(range(1,23)) + ['X', 'Y', 'MT']

CHROMOSOMES = getChr()

##Target rules

rule svs:
    input:
        expand(f"{OUTDIR}/qc/coverage/{{sample}}.stats",
               sample=config["samples"]),
        expand(f"{OUTDIR}/qc/stats/{{sample}}/{{sample}}.stats",
               sample=config["samples"]),
        expand(f"{OUTDIR}/qc/read_length/{{sample}}.txt",
               sample=config["samples"]),
        f"{OUTDIR}/all_annotated/svs_vcfanno.ovl.tab",

rule snps:
    input:
        expand(f"{OUTDIR}/longshot_vep_annotated/all-{{chromosome}}.snps.annot.vcf.gz",
            chromosome=CHROMOSOMES)
