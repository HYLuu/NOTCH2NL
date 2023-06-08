#!/bin/bash

set -e
set -u
set -o pipefail

# Sort all the .proc files for multijoin

for i in *.proc
	do
	base=$(basename ${i} .proc)
	sort ${base}.proc > ${base}.sort
	done
