#!/bin/bash

set -e
set -u
set -o pipefail

for i in *.gz
	do
	bgzip -d $i
	done
