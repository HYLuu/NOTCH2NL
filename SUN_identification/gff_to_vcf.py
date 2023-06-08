'''
This program turns the ‹prefix›.gff format files back into ‹prefix›.vcf format files.
Other than the ‹prefix›.gff format file, the program also needs the reference sequence to be
in the same directory, here is hg38 NOTCH2. Because .gff and .vcf files use slightly different ways
in keeping reference  and query bases.

Input:
‹prefix›.gff - .gff file to convert to .vcf file
Output:
‹prefix›_transform.vcf - .vcf format file of the input ‹‹prefix›.gff file
'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

data = []
final_snp = []
final_indel = []

# Read .gff file
with open(args.filename, 'r') as f:
    for line in f.readlines():
        data.append(line.rstrip('\n').split('\t'))

# Read reference
fasta = []
with open('NOTCH2_hg38_gene.fasta', 'r') as f:
    for line in f.readlines():
        fasta.append(line.rstrip('\n'))
# Merge sequence
seq = "".join([str(item) for item in fasta])

del data[:2]

# Split data and keep only relevant items
for item in data:
    add = item[8].split(';')
    del item[8]
    item.extend(add)

for item in data:
    del item[:3]

# Remove and reorder
indices = [2, 3, 4, 5, 7, 8, 9, 13]
for item in data:
    for i in sorted(indices, reverse=True):
        del item[i]

order = [0, 1, 2, 5, 4, 3]
for item in data:
    item[:] = [item[i] for i in order]

# Convert each line of variant to .vcf format based on the variant type
idx = [2, 3, 6]
for item in data:
    if item[2] == 'Name=substitution':
        item[3] = item[3][10:]
        item[4] = item[4][12:]
        item.insert(0, "hg38_dna_NOTCH2")
        for i in sorted(idx, reverse=True):
            del item[i]
        item.insert(2, ".")
        item.insert(5, ".")
        item.insert(6, "PASS")
        item.insert(7, ".")
        item.insert(8, ".")
        item.insert(9, ".")
        final_snp.append(item)
    if item[2] == 'Name=insertion':
        item[4] = seq[int(item[0]) - 1] + item[4][12:]
        item[3] = seq[int(item[0]) - 1]
        item.insert(0, "hg38_dna_NOTCH2")
        for i in sorted(idx, reverse=True):
            del item[i]
        item.insert(2, ".")
        item.insert(5, ".")
        item.insert(6, "PASS")
        item.insert(7, ".")
        item.insert(8, ".")
        item.insert(9, ".")
        final_indel.append(item)
    if item[2] == 'Name=deletion':
        item[4] = seq[int(item[0]) - 2]
        item[3] = seq[int(item[0]) - 2] + item[3][10:]
        item[0] = str(int(item[0]) - 1)
        item.insert(0, "hg38_dna_NOTCH2")
        for i in sorted(idx, reverse=True):
            del item[i]
        item.insert(2, ".")
        item.insert(5, ".")
        item.insert(6, "PASS")
        item.insert(7, ".")
        item.insert(8, ".")
        item.insert(9, ".")
        final_indel.append(item)

# Output file
with open(args.filename.rstrip('.gff') + '_transform.vcf', 'w') as myfile:
    myfile.write('##fileformat=VCFv4.2' + '\n')
    myfile.write('##source=NucDiffv2.0' + '\n')
    myfile.write('#CHROM' + '\t' + 'POS' + '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + 'FILTER' +  '\t' + 'INFO' + '\t' + 'FORMAT' + '\t' + '${sm}' + '\n')
    for item in data:
        myfile.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\t' + item[6] + '\t' + item[7] + '\t' + item[8] + '\t' + item[9] + '\n')