#!/bin/bash

set -e
set -u
set -o pipefail

for i in *.fasta
	do
	base=$(basename ${i} .fasta)
	# output will be in the directory 'nuc_var/results'
	nucdiff --vcf yes NOTCH2_hg38_gene ${i} nuc_var ${i}_hg38_NOTCH2
	done

