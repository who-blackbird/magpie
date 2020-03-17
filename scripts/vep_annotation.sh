#!/usr/bin/env bash
##Author: Alba Sanchis, as2635@cam.ac.uk

IN=$1
OUT=$2
CHROM=$3

VEPDIR=/rds/project/who1000/rds-who1000-wgs10k/WGS10K/data/projects/nanopore/us/resources/vep
CADDSNPS=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources/CADD/GRCh38/v1.5/whole_genome_SNVs.tsv.gz
CADDINDELS=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources/CADD/GRCh38/v1.5/InDels.tsv.gz
EXACPLI=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources/gnomad/gnomad.v2.1.1.lof_metrics.by_transcript_forVEP.txt
REVELPLUG=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources/REVEL/grch38/new_tabbed_revel.sorted.tsv.gz
GNOMAD38R3=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources//gnomad/GRCh38/r3.0/byChr/gnomad.genomes.r3.0.sites.chr${CHROM}.vcf.bgz
GNOMAD38R2G=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources//gnomad/GRCh38/r2.1.1/byChr/gnomad.genomes.r2.1.sites.grch38.chr${CHROM}_noVEP.vcf.gz
GNOMAD38R2E=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources//gnomad/GRCh38/r2.1.1/byChr/gnomad.exomes.r2.1.sites.grch38.chr${CHROM}_noVEP.vcf.gz
TOPMED=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources//topmed/GRC38/freeze5/bravo-dbsnp-all_noChr.vcf.gz
CLINVAR=/rds/project/flr24/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/resources//clinvar/hg38/20191202/clinvar_20191202.vcf.gz

OUTTMP=${OUT%vcf.gz}TMP.vcf.gz

#VEP
echo "Running VEP"
vep --vcf -i $IN -o STDOUT --format vcf --offline --cache --dir_cache $VEPDIR --force_overwrite \
	--species homo_sapiens \
	--assembly GRCh38 \
	--vcf_info_field ANN \
	--sift b --polyphen b --humdiv --regulatory --allele_number --total_length --numbers --domains \
    --hgvs --protein --symbol --ccds --uniprot --tsl --canonical --mane --biotype --check_existing \
    --af --af_1kg --pubmed --gene_phenotype --variant_class \
	--plugin CADD,${CADDSNPS},${CADDINDELS} \
	--plugin ExACpLI,${EXACPLI} \
	--plugin REVEL,${REVELPLUG} \
	--plugin SpliceRegion \
	--stats_text \
	| bgzip -c > ${OUT}
	
tabix -p vcf $OUT

#Gnomad b38 v3
echo "Annotating with gnomAD v3"

bcftools annotate -c INFO/GNOMADgV3_AC:=INFO/AC,INFO/GNOMADgV3_AF:=INFO/AF,INFO/GNOMADgV3_nhomalt:=INFO/nhomalt,INFO/GNOMADgV3_AC_male:=INFO/AC_male,INFO/GNOMADgV3_AC_female:=INFO/AC_female,INFO/GNOMADgV3_AF_male:=INFO/AF_male,INFO/GNOMADgV3_AF_female:=INFO/AF_female,INFO/GNOMADgV3_AF_afr:=INFO/AF_afr,INFO/GNOMADgV3_AF_amr:=INFO/AF_amr,INFO/GNOMADgV3_AF_asj:=INFO/AF_asj,INFO/GNOMADgV3_AF_eas:=INFO/AF_eas,INFO/GNOMADgV3_AF_fin:=INFO/AF_fin,INFO/GNOMADgV3_AF_nfe:=INFO/AF_nfe,INFO/GNOMADgV3_AF_oth:=INFO/AF_oth,INFO/GNOMADgV3_AF_sas:=INFO/AF_sas,INFO/GNOMADgV3_AF_ami:=INFO/AF_ami,INFO/GNOMADgV3_AC_afr:=INFO/AC_afr,INFO/GNOMADgV3_AC_amr:=INFO/AC_amr,INFO/GNOMADgV3_AC_asj:=INFO/AC_asj,INFO/GNOMADgV3_AC_eas:=INFO/AC_eas,INFO/GNOMADgV3_AC_fin:=INFO/AC_fin,INFO/GNOMADgV3_AC_nfe:=INFO/AC_nfe,INFO/GNOMADgV3_AC_oth:=INFO/AC_oth,INFO/GNOMADgV3_AC_sas:=INFO/AC_sas,INFO/GNOMADgV3_AC_ami:=INFO/AC_ami,INFO/GNOMADgV3_nhomalt_afr:=INFO/nhomalt_afr,INFO/GNOMADgV3_nhomalt_amr:=INFO/nhomalt_amr,INFO/GNOMADgV3_nhomalt_asj:=INFO/nhomalt_asj,INFO/GNOMADgV3_nhomalt_eas:=INFO/nhomalt_eas,INFO/GNOMADgV3_nhomalt_fin:=INFO/nhomalt_fin,INFO/GNOMADgV3_nhomalt_nfe:=INFO/nhomalt_nfe,INFO/GNOMADgV3_nhomalt_oth:=INFO/nhomalt_oth,INFO/GNOMADgV3_nhomalt_sas:=INFO/nhomalt_sas,INFO/GNOMADgV3_nhomalt_ami:=INFO/nhomalt_ami \
	-a $GNOMAD38R3 $OUT -Oz -o $OUTTMP

mv $OUTTMP $OUT

tabix -p vcf $OUT

#Gnomad b38 v2 genomes
echo "Annotating with gnomAD v2 genomes"

bcftools annotate -c INFO/GNOMADg_AC:=INFO/AC,INFO/GNOMADg_AF:=INFO/AF,INFO/GNOMADg_nhomalt:=INFO/nhomalt,INFO/GNOMADg_AC_male:=INFO/AC_male,INFO/GNOMADg_AC_female:=INFO/AC_female,INFO/GNOMADg_AF_male:=INFO/AF_male,INFO/GNOMADg_AF_female:=INFO/AF_female,INFO/GNOMADg_AF_afr:=INFO/AF_afr,INFO/GNOMADg_AF_amr:=INFO/AF_amr,INFO/GNOMADg_AF_asj:=INFO/AF_asj,INFO/GNOMADg_AF_eas:=INFO/AF_eas,INFO/GNOMADg_AF_fin:=INFO/AF_fin,INFO/GNOMADg_AF_nfe:=INFO/AF_nfe,INFO/GNOMADg_AF_oth:=INFO/AF_oth,INFO/GNOMADg_AC_afr:=INFO/AC_afr,INFO/GNOMADg_AC_amr:=INFO/AC_amr,INFO/GNOMADg_AC_asj:=INFO/AC_asj,INFO/GNOMADg_AC_eas:=INFO/AC_eas,INFO/GNOMADg_AC_fin:=INFO/AC_fin,INFO/GNOMADg_AC_nfe:=INFO/AC_nfe,INFO/GNOMADg_AC_oth:=INFO/AC_oth,INFO/GNOMADg_nhomalt_afr:=INFO/nhomalt_afr,INFO/GNOMADg_nhomalt_amr:=INFO/nhomalt_amr,INFO/GNOMADg_nhomalt_asj:=INFO/nhomalt_asj,INFO/GNOMADg_nhomalt_eas:=INFO/nhomalt_eas,INFO/GNOMADg_nhomalt_fin:=INFO/nhomalt_fin,INFO/GNOMADg_nhomalt_nfe:=INFO/nhomalt_nfe,INFO/GNOMADg_nhomalt_oth:=INFO/nhomalt_oth,INFO/GNOMADg_popmax:=INFO/popmax,INFO/GNOMADg_AC_popmax:=INFO/AC_popmax,INFO/GNOMADg_AN_popmax:=INFO/AN_popmax,INFO/GNOMADg_AF_popmax:=INFO/AF_popmax \
	-a $GNOMAD38R2G $OUT -Oz -o $OUTTMP

mv $OUTTMP $OUT

tabix -p vcf $OUT

#Gnomad b38 v2 genomes
echo "Annotating with gnomAD v2 exomes"

bcftools annotate -c INFO/GNOMADe_AC:=INFO/AC,INFO/GNOMADe_AF:=INFO/AF,INFO/GNOMADe_nhomalt:=INFO/nhomalt,INFO/GNOMADe_AC_male:=INFO/AC_male,INFO/GNOMADe_AC_female:=INFO/AC_female,INFO/GNOMADe_AF_male:=INFO/AF_male,INFO/GNOMADe_AF_female:=INFO/AF_female,INFO/GNOMADe_AF_afr:=INFO/AF_afr,INFO/GNOMADe_AF_amr:=INFO/AF_amr,INFO/GNOMADe_AF_asj:=INFO/AF_asj,INFO/GNOMADe_AF_eas:=INFO/AF_eas,INFO/GNOMADe_AF_fin:=INFO/AF_fin,INFO/GNOMADe_AF_nfe:=INFO/AF_nfe,INFO/GNOMADe_AF_oth:=INFO/AF_oth,INFO/GNOMADe_AF_sas:=INFO/AF_sas,INFO/GNOMADe_AC_afr:=INFO/AC_afr,INFO/GNOMADe_AC_amr:=INFO/AC_amr,INFO/GNOMADe_AC_asj:=INFO/AC_asj,INFO/GNOMADe_AC_eas:=INFO/AC_eas,INFO/GNOMADe_AC_fin:=INFO/AC_fin,INFO/GNOMADe_AC_nfe:=INFO/AC_nfe,INFO/GNOMADe_AC_oth:=INFO/AC_oth,INFO/GNOMADe_AC_sas:=INFO/AC_sas,INFO/GNOMADe_nhomalt_afr:=INFO/nhomalt_afr,INFO/GNOMADe_nhomalt_amr:=INFO/nhomalt_amr,INFO/GNOMADe_nhomalt_asj:=INFO/nhomalt_asj,INFO/GNOMADe_nhomalt_eas:=INFO/nhomalt_eas,INFO/GNOMADe_nhomalt_fin:=INFO/nhomalt_fin,INFO/GNOMADe_nhomalt_nfe:=INFO/nhomalt_nfe,INFO/GNOMADe_nhomalt_oth:=INFO/nhomalt_oth,INFO/GNOMADe_nhomalt_sas:=INFO/nhomalt_sas,INFO/GNOMADe_popmax:=INFO/popmax,INFO/GNOMADe_AC_popmax:=INFO/AC_popmax,INFO/GNOMADe_AN_popmax:=INFO/AN_popmax,INFO/GNOMADe_AF_popmax:=INFO/AF_popmax \
	-a $GNOMAD38R2E $OUT -Oz -o $OUTTMP

mv $OUTTMP $OUT

tabix -p vcf $OUT

#TOPMED
echo "Annotating with Topmed"

bcftools annotate -c INFO/TOPMED_AC:=INFO/AC,INFO/TOPMED_AF:=INFO/AF,INFO/TOPMED_AN:=INFO/AN,INFO/TOPMED_Het:=INFO/Het,INFO/TOPMED_Hom:=INFO/Hom \
	-a $TOPMED $OUT -Oz -o $OUTTMP

mv $OUTTMP $OUT

tabix -p vcf $OUT

#CLINVAR
echo "Annotating with Clinvar"

bcftools annotate -c INFO \
	-a $CLINVAR $OUT -Oz -o $OUTTMP

mv $OUTTMP $OUT

tabix -p vcf $OUT
