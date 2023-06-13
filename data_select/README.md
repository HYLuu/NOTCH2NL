# Data Selection

## Data

Genbank data from [HPP_Year1_Assemblies](https://github.com/human-pangenomics/HPP_Year1_Assemblies/tree/main/genbank_changes): <br />

- **contig_translate** is a table of sample contig names and Genbank contig names. Used in `process_NOTCH2_seq.py` and `process_NOTCH2NL_seq.py`. <br />
- **sample_accession** is a table of sample IDs and Genbank accession IDs. Used in `process_NOTCH2_seq.py` and `process_NOTCH2NL_seq.py`. <br />

Extracted data from [Ensembl gene annotations](https://projects.ensembl.org/hprc/) (.gff): <br />

- **NOTCH2NL_annotation** is a subset of NOTCH2NL gene annotations from all the HPRC assembly annotations (.gff3). Used in `process_NOTCH2NL_seq.py`. <br />
- **NOTCH2_R_anno** is a subset of NOTCH2 and NOTCH2NLR genes annotations from all the HPRC assembly annotations (.gff3). Used in `process_NOTCH2_seq.py`. <br />

Sample list: <br />

- **target_assembly** is a sample list of the input of `process_NOTCH2NL_seq.py` and `process_NOTCH2_seq.py`.

## Python scripts

### `Check_NOTCH2_R.py`

```
$ python3 Check_NOTCH2_R.py
```
selects a set of HPRC assemblies that have NOTCH2 and NOTCH2NLR on the same contig but not overlapped. The command outputs the assembly IDs that NOTCH2 and NOTCH2NLR not on the same contig.

### `process_NOTCH2NL_seq.py`

```
$ python3 process_NOTCH2NL_seq.py [list of assembly IDs]
```
produces BED files and a bash script. The bash script is used to extract contig with NOTCH2NLA/B/C.
The bed files are supposed to be used with `bedtools getfasta` for subsetting each NOTCH2NL sequence from the retrieved contig.

### `process_NOTCH2_seq.py`

```
$ python3 process_NOTCH2_seq.py [list of assembly IDs]
```
produces BED files and bash script. The bash script is used to extract contig with NOTCH2 and NOTCH2NLR.
The bed files are supposed to be used with `bedtools getfasta` for subsetting each NOTCH2/NOTCH2NLR sequence from the retrieved contig.

## Bash scripts

### `search_NL_from_gff.sh`

```
$ bash search_NL_from_gff.sh 
```
extracts all the NOTCH2NL gene annotations from all the Ensembl gene annotations and outputs **NOTCH2NL_annotation**.

### `fasta_subset_NL_from_bed.sh`

```
$ bash fasta_subset_NL_from_bed.sh
```
used `bedtools getfasta` and BED files produed from `process_NOTCH2_seq.py` and `process_NOTCH2NL_seq.py` to extract NOTCH2NL sequences from each assembly contig. Outputs <prefix>_ABC.fa files.

### `splitFasta.sh`

```
$ bash splitFasta.sh [file to split]
```
splits sequences in the input fasta file into files. Each output file contains only one sequence.

### `split_ABC.sh`

```
$ bash split_ABC.sh
```
applies `splitFasta.sh` on all <prefix>_ABC.fa files.

## Workflow

[Workflow of data selection step](workflow.png)

![alt text](https://github.com/HYLuu/NOTCH2NL/blob/main/data_select/workflow.png)
