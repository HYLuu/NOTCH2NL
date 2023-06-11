#!/bin/bash

set -e
set -u
set -o pipefail


for i in *_isec; do
	if [ -d "$i" ]; then
		cd "$i"
		base=$(basename ${i} _isec)
		searchstring1="only_"
		head=$(echo "${base#*$searchstring1}" | head -c2)
		if [ "$head" == "HG" ]; then
			searchstring2="HG"
			prefix=${base%$searchstring2*}
			python3 ../../../../hg38_NOTCH2/ratio_sun/clear_multialleic.py ${i}.table ../${prefix}transform.vcf
			sed '1d' ${i}_filtered.table | awk -F"," '$1=$1' OFS="\t" | awk '{$8=($6/$7)}1' OFS="\t" > ${base}_calculated.table
		elif [ "$head" == "NA" ]; then
			searchstring2="NA"
			prefix=${base%$searchstring2*}
			python3 ../../../../hg38_NOTCH2/ratio_sun/clear_multialleic.py ${i}.table ../${prefix}transform.vcf
			sed '1d' ${i}_filtered.table | awk -F"," '$1=$1' OFS="\t" | awk '{$8=($6/$7)}1' OFS="\t" > ${base}_calculated.table
		fi
		cd ..
	fi	
done


