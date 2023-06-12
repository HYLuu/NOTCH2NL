# Copy Number Estimation

## Data

- **SUN_ratio_plot** contains the SUN fraction ratio plots produced in this project. <br />
- **NOTCH2_hg38_gene.fasta** is the reference fasta file used for SUN identification and copy number estimation in this sproject. <br />

## Python scripts

### `clear_multialleic.py`

```
$ python3 clear_multialleic.py [<prefix>.table] [<NOTCH2NL_SUN>_transform.vcf]
```
clears the multialleic variants in the table produced by `GATK VariantsToTable`.

### `plot_ratio.py`

```
$ python3 plot_ratio.py [SUN fraction list of NOTCH2NLA/B] [SUN fraction list of NOTCH2NLC] [SUN fraction list of NOTCH2NLR] [SUN fraction list of unique NOTCH2NLR] [SUN fraction list of NOTCH2]
```
plots the SUN ratio of each individual. Read SUN fractions of each SUN set in th order: NOTCH2NLA/B, NOTCH2NLC, NOTCH2NLR, unique NOTCH2NLR, NOTCH2. Samnple output plots are in directory **SUN_ratio_plot**.

## Bash scripts

### `extract_split_align.sh`

```
$ bash extract_split_align.sh
```
extracts reads aligned to NOTCH2 and NOTCH2NL regions from the CRAM files, splits the reads based on read groups, converts the reads into FASTQ, and aligns them to GRCh38 NOTCH2 using BWA-MEM2.

### `merge_call_variants.sh`

```
$ bash merge_call_variants.sh
```
merges the aligned BAM produced by `extract_split_align.sh` and calls variants from the merged BAMs.

### `vcf_isec.sh`

```
$ bash vcf_isec.sh
```
uses `bcftools isec` for intersecting final VCF produced by `merge_call_variants.sh` and SUN sets VCF of NOTCH2/NOTCH2NL.

### `VariantsToTable.sh`

```
$ bash VariantsToTable.sh
``` 
extracts information needed to calculate SUN fractions in every intersection VCF.

### `fil_calculate_ratio.sh`

```
$ bash fil_calculate_ratio.sh
```
goes into each intersection directory, clears multialleic allels only keeping the correct ones using `clear_multialleic.py`, and output calculated SUN ratio.

### `data_to_plot.sh`

```
$ bash data_to_plot.sh
```
extracts information needed for plotting and merge the SNP and INDEL of the same individual. Output of the script is used to make SUN ratio plots.

## Workflow
