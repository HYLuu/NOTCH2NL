#!/bin/bash

set -e
set -u
set -o pipefail

for i in all_A_B_only*_pos_all.data
	do
	searchstring="only_"
	sample=$(echo "${i#*$searchstring}" | head -c7)
	python3 plot_ratiio.py all_A_B_only_${sample}_pos_all.data clustered_NLC_only_${sample}_pos_all.data NLR_48_only_${sample}_pos_all.data uni6_R_only_${sample}_pos_all.data NL2_all_only_${sample}_pos_all.data
	done
	
