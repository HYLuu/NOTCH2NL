#!/bin/bash

set -e
set -u
set -o pipefail

bash get_new_seq_NLR.sh
rm *.f1_assembly_v2.fa *.f1_assembly_v2.fa.fai
bash fasta_subset_from_bed.sh
