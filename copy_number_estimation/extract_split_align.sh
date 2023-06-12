#!/bin/bash

set -e
set -u
set -o pipefail

# Extract and index reads mapped to NOTCH2 and NOTCH2NL regions (based on GRCh38 on UCSC genome browser) from the CRAM files derived from 1000GP data portal (https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).
for i in *.final.cram
	do
	base=$(basename ${i} .final.cram)
	samtools index $i
	samtools view -O bam -o ${base}_N2.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa $i "chr1:119911553-120069662"
	samtools view -O bam -o ${base}_A.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa $i "chr1:146146521-146229041"
	samtools view -O bam -o ${base}_B.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa $i "chr1:148600285-148679746"
	samtools view -O bam -o ${base}_C.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa $i "chr1:149390621-149471833"
	samtools view -O bam -o ${base}_R.bam -T GRCh38_full_analysis_set_plus_decoy_hla.fa $i "chr1:120723912-120794851"
	samtools index ${base}_N2.bam
	samtools index ${base}_A.bam
	samtools index ${base}_B.bam
	samtools index ${base}_C.bam
	samtools index ${base}_R.bam
	done

# For each BAM file gotten from last step, split them into several BAM files based on read groups. This step is needed for marking duplication in GATK.
for j in *.bam
	do
	base=$(basename ${j} .bam)
	# Get read group ID for each BAM
	samtools view -H ${j} | grep '^@RG' | awk '{print substr($2,4), substr($9,4), substr($7,4), substr($3,4), substr($5,4)}' > ${base}.RGID
	IN_FILE="./${base}.RGID"
	# Split based on read group
	while read -ra LINE
		do
		a="${LINE[0]}"
		samtools view -b -r $a ${j} > ${base}_${a}_spt.bam
		done < "$IN_FILE"
	done

# Convert the split BAM files into FASTQ files for realignment.
for k in *_spt.bam
	do
	base=$(basename ${k} _spt.bam)
	# Sorting
	samtools sort -n -o ${base}.qsort.bam ${k}
	# BAM to FASTQ
	bedtools bamtofastq -i ${base}.qsort.bam -fq ${base}_1.fq -fq2 ${base}_2.fq
	done

# Map FASTQ to reference GRCh38 NOTCH2 and provide read group information. This information is key for downstream GATK functionality. GATK will not work without a read group tag.
for m in *.qsort.bam
	do
	base=$(basename ${m} .qsort.bam)
	ID=$(echo $base | cut -d_ -f3-6)
	RG=$(echo $base | cut -d_ -f1-2)
	line=$(awk '$1 ~ /^${ID}$/' ${RG}.RGID)
	read id pu sm pl lb <<< ${ID}
	~/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -K 100000000 -Y -R '@RG\tID:${id}\tPU:${pu}\tSM:${sm}\tPL:${pl}\tLB:${lb}' NOTCH2_hg38_gene.fasta ${base}_1.fq ${base}_2.fq > ${base}_hg38_NOTCH2.sam
	done
	
