#!/bin/bash

set -e
set -u
set -o pipefail

for i in *.gff3
	do
	base=$(basename ${i} .gff3)
	awk 'BEGIN{print ARGV[1]} /NOTCH2NL/ {print $1,$4,$5,$7,$9}' ${i} >> NOTCH2NL_annotation
	done
