#!/bin/bash

set -e
set -u
set -o pipefail

#mkdir process_to_plot/

for i in *_isec; do
	if [ -d "$i" ]; then
		cd "$i"
		base=$(basename ${i} _isec)
		awk '{print $2,$8}' ${base}_calculated.table > ${base}_pos.data
		cp ${base}_pos.data ../process_to_plot/
		cd ..
	fi	
done

cd process_to_plot/

for j in *_filtered_indels_pos.data; do
	base=$(basename ${j} _filtered_indels_pos.data)	
	cat ${base}_filtered_indels_pos.data ${base}_filtered_snps_pos.data > ${base}_pos_all.data
done
