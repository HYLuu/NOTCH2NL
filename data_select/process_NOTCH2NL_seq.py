'''
This program produces BED files and bash script. The bash script is used to extract contig with NOTCH2NLA/B/C.
The bed files are supposed to be used with bedtools getfasta for subsetting each NOTCH2NL sequence from the retrieved contig.

Assembly data from HPRC are used for accession number and contig name conversion.
Data used:
<sample_accession> - table of sample IDs and Genbank accession IDs
<contig_translate> - table of sample contig names and Genbank contig names
(Above are from: https://github.com/human-pangenomics/HPP_Year1_Assemblies/tree/main/genbank_changes)
<NOTCH2NL_annotation> - subset of NOTCH2NL gene annotations from all the HPRC assembly annotations (.gff3)

Input:
List of assemblies to extract NOTCH2NL sequences from. e.g. HG00438#1
Output:
<prefix>.bed - BED files subsetting each NOTCH2NL sequence from contig. Contains NOTCH2NLA/B/C
get_NL_seq_from_assembly.sh - a bash script using samtools faidx to extract contigs from original assemblies
'''

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

id = []
anno = []
contig = []
trans = []

### Read files and convert into dataframes
with open('sample_accession', "r") as myfile:
    for line in myfile.readlines():
        id.append(line.rstrip('\n').split('\t'))
sample_id = pd.DataFrame(id)
with open('NOTCH2NL_annotation', "r") as myfile:
    for line in myfile.readlines():
        anno.append(line.rstrip('.gff3\n').lstrip('homo_sapiens_gca').split(' '))
annotation = pd.DataFrame(anno)
with open('contig_translate', "r") as myfile:
    for line in myfile.readlines():
        trans.append(line.rstrip('\n').split('\t'))
translate = pd.DataFrame(trans)
# Target assemblies
with open(args.filename, 'r') as f:
    for line in f.readlines():
        contig.append(line.rstrip('\n').split('#'))
contig_frame = pd.DataFrame(contig)

# Set index to sample name
contig_frame = contig_frame.set_index(0)
sample_id = sample_id.set_index(0)

### Add a column of Genbank IDs to contig_frame
# Find the corresponding Genbank IDs from sample_id
contig_frame[3] = np.where(contig_frame[1] == '1', sample_id.loc[contig_frame.index][2], sample_id.loc[contig_frame.index][1])
contig_frame[3] = contig_frame[3].str.replace(r'\..*', '', regex=True)
contig_frame[3] = contig_frame[3].str.replace('GCA_', '', regex=True)

# Add a column of 1s into annotation for the bed files
annotation.insert(5, 5, 1)
# Reorder
annotation = annotation[[0, 1, 2, 4, 5, 3]]
# Clear data only keep the necessary information
annotation[0] = annotation[0].str.replace(r'\..*', '', regex=True)
annotation[0] = annotation[0].str.replace(r'v.*', '', regex=True)
annotation[4] = annotation[4].str.replace(r';biotype=protein_coding;description=notch', '', regex=True)
annotation[4] = annotation[4].str.replace(r'.*.\;', '', regex=True)
annotation[4] = annotation[4].str.replace('Name=', '', regex=True)
# Replace the Genbank contig names with sample contig names
annotation[0] = annotation[0].map(lambda x: translate.set_index(1)[0].get(x, x))
# Remove the rows of NOTCH2NLR
indexR = annotation[annotation[4] == 'NOTCH2NLR'].index
annotation.drop(indexR, inplace=True)

# Set index to Genbank IDs
contig_frame = contig_frame.set_index(3)
# Replace "1" and "2" with "paternal" and "maternal"
contig_frame.loc[contig_frame[1] == "1", [1]] = "paternal"
contig_frame.loc[contig_frame[1] == "2", [1]] = "maternal"

### Output BED files
# Iterate through each target sample in contig_frame
for index, row in contig_frame.iterrows():
    # Get the index of sample in annotation data frame
    indexAnno = annotation[annotation[0] == index].index
    # Get integer location of the sample in annotation data frame because index =/= integer location
    indexAnno_i = annotation.index.get_loc(indexAnno[0])
    # Use sample ID and contig name as BED file name
    with open(str(annotation.iloc[indexAnno_i+1, 0] + ".bed"), "w") as f:
        f.write(str(annotation.iloc[indexAnno_i + 1, 0]) + '\t' + str(annotation.iloc[indexAnno_i + 1, 1]) + '\t' + str(annotation.iloc[indexAnno_i + 1, 2]) + '\t' + str(annotation.iloc[indexAnno_i + 1, 3]) + '\t' + str(annotation.iloc[indexAnno_i + 1, 4]) + '\t' + str(annotation.iloc[indexAnno_i + 1, 5]) + '\n')
        f.write(str(annotation.iloc[indexAnno_i + 2, 0]) + '\t' + str(annotation.iloc[indexAnno_i + 2, 1]) + '\t' + str(annotation.iloc[indexAnno_i + 2, 2]) + '\t' + str(annotation.iloc[indexAnno_i + 2, 3]) + '\t' + str(annotation.iloc[indexAnno_i + 2, 4]) + '\t' + str(annotation.iloc[indexAnno_i + 2, 5]) + '\n')
        f.write(str(annotation.iloc[indexAnno_i + 3, 0]) + '\t' + str(annotation.iloc[indexAnno_i + 3, 1]) + '\t' + str(annotation.iloc[indexAnno_i + 3, 2]) + '\t' + str(annotation.iloc[indexAnno_i + 3, 3]) + '\t' + str(annotation.iloc[indexAnno_i + 3, 4]) + '\t' + str(annotation.iloc[indexAnno_i + 3, 5]) + '\n')

### Output bash script
with open("get_NL_seq_from_assembly.sh", "w") as f:
    # Write header
    f.write('#!/bin/bash\n')
    f.write('\n')
    f.write('set -e\n')
    f.write('set -u\n')
    f.write('set -o pipefail\n')
    f.write('\n')
    # Iterate through each target sample in contig_frame
    for index, row in contig_frame.iterrows():
        indexAnno = annotation[annotation[0] == index].index
        indexAnno_i = annotation.index.get_loc(indexAnno[0])
        f.write('samtools faidx ' + str(annotation.iloc[indexAnno_i+1, 0]).split('#')[0] + '.' + str(row[1]) + '.f1_assembly_v2.fa' + ' ' + str(annotation.iloc[indexAnno_i+1, 0]) + ' > ' + str(annotation.iloc[indexAnno_i+1, 0]) + '.fa' + '\n')
