'''
This program clears the multialleic variants in the table produced by GATK VariantsToTable. The '-SMA' option in
GATK VariantsToTable splits multialleic variants into multiple lines with single alternate allele. This program compares
the input table with the NOTCH2NL SUN vcf and remove the alternate alleles that are not part of the SUNs.

Input:
<prefix>.table - a table produced by GATK VariantsToTable, which contains part of the information in the original .vcf file (a intersection vcf of sample variants and SUNs)
<NOTCH2NL_SUN>_transform.vcf - .vcf file of the SUN set used for intersection vcf of the above <prefix>.table
Output:
<prefix>_filtered.table - the original table with multialleic variants cleared
'''


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('reference')
args = parser.parse_args()

# Read <prefix>.table
data = []
with open(args.filename, 'r') as f:
    for line in f.readlines():
        data.append(line.rstrip('\n').split('\t'))

# Read <NOTCH2NL_SUN>_transform.vcf
ref = []
with open(args.reference, 'r') as f:
    for line in f.readlines():
        ref.append(line.rstrip('\n').split('\t'))

del data[:1]
del ref[:3]

i = 0
rf = 0
to_del = []
# Walk through each line of variants
while i < len(data):
    # If the variant is multialleic
    if data[i][4].count(',') > 1:
        for j in range(len(ref)):
            # Find the corresponding variant in <NOTCH2NL_SUN>_transform.vcf and save its index in variable rf
            if data[i][1] == ref[j][1]:
                rf = j
        t = 0  # index of AD
        while data[i][1] == ref[rf][1]:  # For all the alternate alleles in a multialleic variant
            if i == (len(data) - 1):  # In the last two line of data
                if data[i][2] != ref[rf][3] or data[i][3] != ref[rf][4]:  # The allele that is not in the SUN set
                    to_del.append(i)  # Save to to-be-removed
                    t += 1  # Move to next AD
                elif data[i][2] == ref[rf][3] and data[i][3] == ref[rf][4]:  # The allele is in the SUN set
                    data[i][4] = ','.join([data[i][4].split(',')[0], data[i][4].split(',')[t + 1]])  # Remove the AD of the allele that is not in the SUN set
                i += 1  # Move to next line
                break
            elif i < (len(data) - 1):  # Not in the last two line of data
                if data[i][2] != ref[rf][3] or data[i][3] != ref[rf][4]:  # The allele that is not in the SUN set
                    to_del.append(i)  # Save to to-be-removed
                    t += 1  # Move to next AD
                    i += 1  # Move to next line
                elif data[i][2] == ref[rf][3] and data[i][3] == ref[rf][4]:  # The allele is in the SUN set
                    data[i][4] = ','.join([data[i][4].split(',')[0], data[i][4].split(',')[t + 1]])  # Remove the AD of the allele that is not in the SUN set
                    i += 1  # Move to next line
    else:  # If the variant is not multialleic
        i += 1

# Delete the lines in to-be-removed
for index in sorted(to_del, reverse=True):
    del data[index]

# Output the filtered table
with open(args.filename.rstrip('.table') + '_filtered.table', 'w') as myfile:
    myfile.write('#CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + '${sm}.AD' + '\t' + '${sm}.DP' + '\n')
    for item in data:
        myfile.write(
            item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')