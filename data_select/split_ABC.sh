#!/bin/bash

set -e
set -u
set -o pipefail

for i in *_ABC.fa
	do
	bash splitFasta.sh $i
	done

