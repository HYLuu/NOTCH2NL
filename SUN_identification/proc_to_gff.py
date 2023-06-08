'''
This program turns the ‹prefix›.proc format files back into ‹prefix›.gff format files.
Other than the ‹prefix›.proc format file, the program also needs the original ‹prefix›.gff file as input.
If the ‹prefix›.proc format file is a joined/excluded file of multiple ‹prefix›.gff files, the input ‹prefix›.gff file
will be the first processed ‹prefix›.gff file.

Input:
‹prefix›.proc - .proc format file
‹original›.gff - original .gff file or the first processed .gff file
Output:
‹prefix›.gff - .gff format file of the input ‹prefix›.proc file, including only variants in the input ‹prefix›.proc file
'''

import argparse

# Read files
parser = argparse.ArgumentParser()
parser.add_argument('filename')  # .proc file
parser.add_argument('reference')  # .gff file
args = parser.parse_args()

data = []
gff = []
head = []

with open(args.filename, 'r') as f:
    for line in f.readlines():
        data.append(line.rstrip('\n').split(' '))
with open(args.reference, 'r') as f:
    for line in f.readlines():
        gff.append(line.rstrip('\n').split('\t'))

# Split items in the .proc file
for item in data:
    add = item[0].split(';')
    del item[0]
    item.extend(add)

# Get the header of .gff file
for k in range(2):
    head.append(gff[k])
del gff[:2]

# Split items in the .gff file
for item in gff:
    add = item[8].split(';')
    del item[8]
    item.extend(add)

# Compare variants in both files. The corresponding line of .proc file in .gff file will be saved for the output .gff.
final = []
for i in range(len(data)):
    for j in range(len(gff)):
        # If the items in .proc file is the same as the one in .gff file
        if data[i][1] == gff[j][3] and data[i][2] == gff[j][4] and data[i][3] == gff[j][9] and data[i][4] == gff[j][-2] and data[i][5] == gff[j][-3] and data[i][0] == gff[j][-4]:
            final.append(gff[j])

# Merge items
for item in final:
    item[8:] = [';'.join(item[8:])]

# Reorder
final = sorted(final, key=lambda x: (int(x[3]), int(x[4])))

# Add header
final.insert(0, head[1])
final.insert(0, head[0])

# Output the .gff format file
with open(args.filename + '.gff', 'w') as myfile:
    myfile.write(final[0][0] + '\n')
    myfile.write(final[1][0] + '\t' + final[1][1] + '\t' + final[1][2] + '\t' + final[1][3] + '\n')
    for i in range(2, len(final)):
        myfile.write(final[i][0] + '\t' + final[i][1] + '\t' + final[i][2] + '\t' + final[i][3] + '\t' + final[i][4] + '\t' + final[i][5] + '\t' + final[i][6] + '\t' + final[i][7] + '\t' + final[i][8] + '\n')