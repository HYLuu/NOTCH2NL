#!/bin/bash

set -e
set -u
set -o pipefail


for i in *.fa
	do
	base=$(basename ${i} .fa)
	bedtools getfasta -fi ${i} -bed ${base}.bed -s -name > ${base}_ABC.fa
	done


