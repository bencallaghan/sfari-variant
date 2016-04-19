# sfari-variant

Pipeline for prioritising variants for SFARI project genes of interest

## Inputs:

Place all input files in folder ./inputs/GENENAME

1. Fasta file of the gene of interest (canonical transcript)
2. .bed file of the gene of interest (canonical transcript)
3. Variants of interest for gene of interest in .vcf format (chrom,start,stop,ref,alt) - MARV variants, Clinvar, whatever
4. Predictprotein scores for the gene

## Scripts:

## Current Work:

Generalize variant prioritisation pipeline in (GENE_prioritised_variants.R) 
