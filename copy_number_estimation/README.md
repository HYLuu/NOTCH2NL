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

### `nucdiff_HPRC.sh`

```
$ nucdiff_HPRC.sh
```
calls variants of each NOTCH2NL sequence (produced in data selection step) )using `Nucdiff` with GRCh38 NOTCH2 as reference.

### `sort_ref.sh`

```
$ bash sort_ref.sh
```
sorts all the .proc files as .sort files. Th sorted files will be used for finding joint SNPs using `multijoin`.

### `multijoin`

```
$ bash multijoin [all the .sort files for multijoin]
```
joins multiple files to find the set of SNPs that shared by all the input files. Standard of joint SNPs could be find in `Nucdiff_gff_process.py`.

## Workflow


