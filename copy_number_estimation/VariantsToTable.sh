#!/bin/bash

set -e
set -u
set -o pipefail

sm_vcf=(*_filtered_indels_final.vcf *_filtered_snps_final.vcf)

for i in *_isec; do
	if [ -d "$i" ]; then
		cd "$i"
		~/gatk-4.4.0.0/gatk VariantsToTable -V 0000.vcf -F CHROM -F POS -F REF -F ALT -SMA true -GF AD -GF DP -O ${i}.table
		cd ..
	fi	
done


