# Distinguish Between NOTCH2NL Genes Using Gene Specific Variants

## Project Description

Assemblies from the Human Pangenome Reference Consortium (HPRC) were used to retrieve singly-unique nucleotides (SUNs) for each NOTCH2NL allele. These SUNs were used to estimate the NOTCH2NL gene copy number from short-read sequencing data generated for 1000 genome samples included in the pangenome data set. Project details can be found in the [report](https://drive.google.com/file/d/1-UVtCL_IFVTf8fIYiBOZ724giBFTZY3p/view?usp=share_link). The three directories are arranged according to the *Methods* section of the report.

## Directories

The three directories are arranged according to the **Methods** section of the report.

*data_select* contains scripts and data used for HPRC assembly selection and extraction.
*SUN_identification* contains scripts used to get variants from each assembly, process the variants and find the joint SNPs for each NOTCH2NL.
*copy_number_estimation* contains scripts used to get variants from short-read sequencing data, find the intersection of these variants and SUNs, calculate SUN fractions, and produce SUN fraction plots. This directory also includes SUN ratio plots from 39 short-read sequencing data samples.


