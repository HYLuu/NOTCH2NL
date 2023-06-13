# SUN Identification

## Python scripts

### `Nucdiff_gff_process.py`

```
$ python3 Nucdiff_gff_process.py [GFF file]
```
turns the ‹prefix›_ref_snps.gff files produced by Nucdiff to ‹prefix›.proc files. ‹prefix›.proc files are used to find the SUNs because they only contain necessary simplified information
for comparing to variants of other assemblies.

### `proc_to_gff.py`

```
$ python3 proc_to_gff.py [.proc file] [original GFF file]
```
turns the ‹prefix›.proc format files back into ‹prefix›.gff format files. Other than the ‹prefix›.proc format file, the program also needs the original ‹prefix›.gff file as input. Details of original ‹prefix›.gff file can be found in program comments.

### `gff_to_vcf.py`

```
$ python3 gff_to_vcf.py [GFF file]
```
turns the ‹prefix›.gff format files back into ‹prefix›.vcf format files. Other than the ‹prefix›.gff format file, the program also needs the reference sequence to be
in the same directory, here is hg38 NOTCH2.

## Bash scripts

### `nucdiff_HPRC.sh`

```
$ bash nucdiff_HPRC.sh
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

[Workflow of SUN identification step](workflow.png)

![alt text](https://github.com/HYLuu/NOTCH2NL/blob/main/SUN_identification/workflow.png)
