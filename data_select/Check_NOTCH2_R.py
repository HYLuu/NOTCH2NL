'''
This script selects a set of HPRC assemblies that have NOTCH2 and NOTCH2NLR on the same contig but not overlapped.

Input data used:
NOTCH2_R_anno - the annotations of NOTCH2 and NOTH2NLR extracted from all the HPRC assembly annotations
NLR_1st_edition.csv - a list of assemblies used for getting SUNs of NOTCH2NLA, NOTCH2NLB and NOTCH2NLC

Output:
Print the assembly names that NOTCH2 and NOTCH2NLR not on the same contig
1st_edit_NLR - list of assembly numbers of assemblies used for getting SUNs of NOTCH2NLA, NOTCH2NLB and NOTCH2NLC
new_edit_NLR - list of assembly numbers of assemblies used for getting SUNs of NOTCH2 and NOTCH2NLR
'''

import pandas as pd

'''
Function for checking whether the coordinates of two sequences are overlapped.
Input: coordinates of two sequences: (a, b) and (x, y) while b > a and y > x
Output: True if overlapped; False if not.
'''
def check_overlap(a, b, x, y):
    if a <= x < b <= y:
        return True
    elif x <= a < y <= b:
        return True
    elif a <= x < y <= b:
        return True
    elif x <= a < b <= y:
        return True


### Get all NLR and NOTCH2NL and check if overlap ###

anno = []
id = []

# Read the annotations of NOTCH2 and NOTCH2NLR extracted from all the HPRC assemblies
with open('NOTCH2_R_anno', "r") as myfile:
    for line in myfile.readlines():
        anno.append(line.rstrip('.gff3\n').split(' '))
annotation = pd.DataFrame(anno)

# Only remain useful information: assembly name, contig, coordinates, gene
annotation[4] = annotation[4].str.replace(r';biotype=protein_coding;description=notch', '', regex=True)
annotation[4] = annotation[4].str.replace(r'.*.\;', '', regex=True)
annotation[4] = annotation[4].str.replace('Name=', '', regex=True)

# Get the first column (assembly names and contigs)
name = annotation[0]

# Check if the NOTCH2 and NOTCH2NLR in each assembly are on the same contig
for i in range(len(name) - 2):
    if name[i].startswith('Homo') and i == len(name) - 3:
        if name[i + 1] != name[i + 2]:
            print(name[i], " NOTCH2 and NOTCH2NLR not on the same contig.")
    elif name[i].startswith('Homo') and name[i + 3].startswith('Homo'):
        if name[i + 1] != name[i + 2]:
            print(name[i], " NOTCH2 and NOTCH2NLR not on the same contig.")
    elif name[i].startswith('Homo') and not name[i + 3].startswith('Homo'):
        ct = 1
        while not name[i + ct + 1].startswith('Homo'):
            if name[i + ct] != name[i + ct + 1]:
                print(name[i], " NOTCH2 and NOTCH2NLR not on the same contig.")
            ct += 1

'''
Results:
Homo_sapiens-GCA_018469415.1-2022_07-genes  NOTCH2 and NOTCH2NLR not on the same contig.
Homo_sapiens-GCA_018469415.1-2022_07-genes  NOTCH2 and NOTCH2NLR not on the same contig.
Homo_sapiens-GCA_018852605.1-2022_07-genes  NOTCH2 and NOTCH2NLR not on the same contig.
'''

# Remove the assemblies that NOTCH2 and NOTCH2NLR not on the same contig
indexA = annotation[annotation[0] == 'Homo_sapiens-GCA_018469415.1-2022_07-genes'].index
indexB = annotation[annotation[0] == 'Homo_sapiens-GCA_018852605.1-2022_07-genes'].index
annotation.drop(indexA, inplace=True)
annotation.drop(indexA + 1, inplace=True)
annotation.drop(indexA + 2, inplace=True)
annotation.drop(indexA + 3, inplace=True)
annotation.drop(indexB, inplace=True)
annotation.drop(indexB + 1, inplace=True)
annotation.drop(indexB + 2, inplace=True)

annotation.reset_index(drop=True, inplace=True)

# Get each column
col1 = annotation[1]
col2 = annotation[2]
name = annotation[0]

# Check whether the NOTCH2 and NOTCH2NLR in each assembly are overlapped
not_overlap = []
k = 0  # Number of NOTCH/NOTCH2NLR non-overlapped assemblies
for i in range(1, len(col1) - 1, 3):  # For each assembly
    a = int(col1[i])  # Get coordinates
    b = int(col2[i])
    x = int(col1[i + 1])
    y = int(col2[i + 1])
    if not check_overlap(a, b, x, y):
        k += 1
        not_overlap.append(name[i - 1])

NLR_valid = pd.DataFrame(not_overlap)

### Compare to the assembly selection of NOTCH2NLA, NOTCH2NLB and NOTCH2NLC ###

# Read the assembly list of the assembly selection of NOTCH2NLA, NOTCH2NLB and NOTCH2NLC
first_edt = pd.read_csv("NLR_1st_edition.csv")

# Clean the data to remain only the accession numbers
NLR_valid[0] = NLR_valid[0].str.replace('Homo_sapiens-GCA_0', '', regex=True)
NLR_valid[0] = NLR_valid[0].str.replace(r'.1-.*', '', regex=True)
first_edt = first_edt['2']

# Output both list of assemblies
with open("1st_edit_NLR", "w") as f:
    for i in range(len(first_edt)):
        f.write(str(first_edt.iloc[i]) + "\n")

with open("new_edit_NLR", "w") as f:
    for i in range(len(NLR_valid)):
        f.write(str(NLR_valid[0].iloc[i]) + "\n")
