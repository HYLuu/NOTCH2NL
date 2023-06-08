'''
This program turns the ‹prefix›_ref_snps.gff files produced by Nucdiff to ‹prefix›.proc files.
‹prefix›.proc files are used to find the SUNs because they only contain necessary simplified information
for comparing to variants of other assemblies. (Necessary for comparison: reference coordinates, variant type,
reference base and query base. Query coordinates are included in the second column as keys for later conversion back to .gff files)

Input:
‹prefix›_ref_snps.gff - .gff file produced by Nucdiff
Output:
‹prefix›_ref_snps.proc - .proc file for variants comparison
'''

import argparse

# Read file
data = []
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

with open(args.filename, 'r') as f:
    for line in f.readlines():
        data.append(line.rstrip('\n').split('\t'))

del data[:2]

# Split the data into items in a list
for item in data:
    add = item[8].split(';')
    del item[8]
    item.extend(add)

# Remain only necessary information
for item in data:
    del item[:3]

indices = [2, 3, 4, 5, 7, 8, 9, 13]
for item in data:
    for i in sorted(indices, reverse=True):
        del item[i]

# Reorder
order = [0, 1, 2, 5, 4, 3]
for item in data:
    item[:] = [item[i] for i in order]

# Merge necessary items for comparison
for item in data:
    item[:5] = [';'.join(item[:5])]

# Output .proc file
with open(args.filename.rstrip('.gff') + '.proc', 'w') as myfile:
    for item in data:
        myfile.write(item[0] + '\t' + item[1] + '\n')