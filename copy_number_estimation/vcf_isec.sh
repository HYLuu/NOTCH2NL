#!/bin/bash

set -e
set -u
set -o pipefail

sm_vcf=(*_filtered_indels_final.vcf *_filtered_snps_final.vcf)
sm_vcf_gz=(*_filtered_indels_final.vcf.gz *_filtered_snps_final.vcf.gz)

for i in ${sm_vcf[@]}
	do
	bgzip -f -c $i > "${i}.gz"
	bcftools index -f ${i}.gz > ${i}.gz.csi
	done

for i in *_transform.vcf
	do
	bgzip -f -c $i > "${i}.gz"
	bcftools index -f ${i}.gz > ${i}.gz.csi
	done

for j in *_transform.vcf.gz
	do
	base1=$(basename ${j} _transform.vcf.gz)
	for k in ${sm_vcf_gz[@]}
		do
		base2=$(basename ${k} _final.vcf.gz)
		bcftools isec -p ${base1}_${base2}_isec -n=2 -c some -w1 $k $j
		done
	done


