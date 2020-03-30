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

CHROMOSOMES = getChr()

# Target rules #

rule structural_variants:
    input:
        expand("{OUTDIR}/minimap2/coverage/{{sample}}.stats",
               sample=config["samples"]),
<<<<<<< HEAD
        f"{OUTDIR}/minimap2/sniffles_annotated/genotypes_annot.vcf.gz",
        f"{OUTDIR}/ngmlr/sniffles_annotated/genotypes_annot.vcf.gz"


rule snps:
    input:
        expand(f"{OUTDIR}/minimap2/longshot_split/{{sample}}-{{chromosome}}.snps.vcf.gz",
            sample=config["samples"], chromosome=CHROMOSOMES)
        #expand(f"{OUTDIR}/minimap2/longshot_vep_annotated/all-{{chromosome}}.snps.annot.vcf.gz",
            #chromosome=CHROMOSOMES)
=======
        expand("{OUTDIR}/minimap2/sniffles_svs/{{sample}}.vcf",
               sample=config["samples"]),
        expand("{OUTDIR}/minimap2/sniffles_annotate/{{sample}}_annot.vcf.tsv",
               sample=config["samples"]),
        expand("{OUTDIR}/minimap2/longshot/{{sample}}.merged.vcf.gz",
               sample=config["samples"])
>>>>>>> c2c1e1d9f811e1aff3311ea9b2f304e708844205
