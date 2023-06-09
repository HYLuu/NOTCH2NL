#!/bin/bash

set -e
set -u
set -o pipefail


while read line ; do
	if [ ${line:0:1} == ">" ] ; then
		index=$(echo "$1" |  cut -d "_" -f1)
    		filename=$(echo "${index}_$line" | cut -d ":" -f1 | tr -d ">")
    		touch ./"$filename".fasta
    		echo "$line" >> ./"${filename}".fasta
  	else
    		echo "$line" >> ./"${filename}".fasta
	fi
done < $1
