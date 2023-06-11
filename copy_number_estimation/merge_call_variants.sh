#!/bin/bash

set -e
set -u
set -o pipefail

#for i in *.sam
#	do
#	base=$(basename ${i} .sam)
#	samtools view -S -b ${i} > ${base}.bam
#	samtools sort ${base}.bam -o ${base}.sorted.bam
#	done

#for j in *.sorted.bam
#	do
#	base=$(basename ${j} .sorted.bam)
#	prefix=$(echo $base | cut -d_ -f1)
#	samtools merge -c -f ${prefix}_all.bam ${prefix}_*.sorted.bam
#	done

#mkdir sep_bam
#mv *_hg38_NOTCH2.sam *_hg38_NOTCH2.sorted.bam *_hg38_NOTCH2.bam sep_bam/

#for k in *_all.bam
#	do
#	base=$(basename ${k} _all.bam)
#	~/gatk-4.4.0.0/gatk MarkDuplicatesSpark -I $k -M ${base}_dedup_metrics.txt -O ${base}_sorted_dedup.bam
#	done

#for m in *_sorted_dedup.bam
#	do
#	base=$(basename ${m} _sorted_dedup.bam)
#	java -jar ~/picard.jar CollectAlignmentSummaryMetrics R=../../../../NOTCH2_hg38_gene.fasta I=$m O=${base}_aln_metrics.txt
#	samtools depth -a $m > ${base}_depth_out.txt
#	~/gatk-4.4.0.0/gatk HaplotypeCaller -R ../../../../NOTCH2_hg38_gene.fasta -I $m -O ${base}_raw_variants.vcf
#	done

#for n in *_raw_variants.vcf
#	do
#	base=$(basename ${n} _raw_variants.vcf)
#	~/gatk-4.4.0.0/gatk SelectVariants -R ../../../../NOTCH2_hg38_gene.fasta -V $n -select-type SNP -O ${base}_raw_snps.vcf
#	~/gatk-4.4.0.0/gatk SelectVariants -R ../../../../NOTCH2_hg38_gene.fasta -V $n -select-type INDEL -O ${base}_raw_indels.vcf
#	~/gatk-4.4.0.0/gatk VariantFiltration -R ../../../../NOTCH2_hg38_gene.fasta -V ${base}_raw_snps.vcf -O ${base}_filtered_snps.vcf \
#	-filter-name "QD_filter" -filter "QD < 2.0" \
#	-filter-name "FS_filter" -filter "FS > 60.0" \
#	-filter-name "MQ_filter" -filter "MQ < 40.0" \
#	-filter-name "SOR_filter" -filter "SOR > 4.0" \
#	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
#	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#	~/gatk-4.4.0.0/gatk VariantFiltration -R ../../../../NOTCH2_hg38_gene.fasta -V ${base}_raw_indels.vcf -O ${base}_filtered_indels.vcf \
#	-filter-name "QD_filter" -filter "QD < 2.0" \
#	-filter-name "FS_filter" -filter "FS > 200.0" \
#	-filter-name "SOR_filter" -filter "SOR > 10.0"
#	~/gatk-4.4.0.0/gatk SelectVariants --exclude-filtered true -V ${base}_filtered_snps.vcf -O ${base}_bqsr_snps.vcf
#	~/gatk-4.4.0.0/gatk SelectVariants --exclude-filtered true -V ${base}_filtered_indels.vcf -O ${base}_bqsr_indels.vcf
#	done

for p in *_sorted_dedup.bam
	do
	base=$(basename ${p} _sorted_dedup.bam)
	~/gatk-4.4.0.0/gatk BaseRecalibrator -R ../../../../NOTCH2_hg38_gene.fasta -I $p --known-sites ${base}_bqsr_snps.vcf --known-sites ${base}_bqsr_indels.vcf -O ${base}_recal_data.table
	~/gatk-4.4.0.0/gatk ApplyBQSR -R ../../../../NOTCH2_hg38_gene.fasta -I $p -bqsr ${base}_recal_data.table -O ${base}_recal_reads.bam
	~/gatk-4.4.0.0/gatk BaseRecalibrator -R ../../../../NOTCH2_hg38_gene.fasta -I ${base}_recal_reads.bam --known-sites ${base}_bqsr_snps.vcf --known-sites ${base}_bqsr_indels.vcf -O ${base}_post_recal_data.table
	~/gatk-4.4.0.0/gatk HaplotypeCaller -R ../../../../NOTCH2_hg38_gene.fasta -I ${base}_recal_reads.bam -O ${base}_raw_variants_recal.vcf
	~/gatk-4.4.0.0/gatk SelectVariants -R ../../../../NOTCH2_hg38_gene.fasta -V ${base}_raw_variants_recal.vcf -select-type SNP -O ${base}_raw_snps_recal.vcf
	~/gatk-4.4.0.0/gatk SelectVariants -R ../../../../NOTCH2_hg38_gene.fasta -V ${base}_raw_variants_recal.vcf -select-type INDEL -O ${base}_raw_indels_recal.vcf
	~/gatk-4.4.0.0/gatk VariantFiltration -R ../../../../NOTCH2_hg38_gene.fasta -V ${base}_raw_snps_recal.vcf -O ${base}_filtered_snps_final.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
	~/gatk-4.4.0.0/gatk VariantFiltration -R ../../../../NOTCH2_hg38_gene.fasta -V ${base}_raw_indels_recal.vcf -O ${base}_filtered_indels_final.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
	done
