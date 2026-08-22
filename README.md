# Ldec_transcriptomics
Reproducible workflows for transcriptomic, functional annotation, immune response, and virus-homology analyses in Leptinotarsa decemlineata.

## Project overview
This repository contains scripts and workflows used to analyse
the transcriptomic response of the Colorado potato beetle
(*Leptinotarsa decemlineata*) to infection with the entomopathogenic fungi
*Beauveria bassiana* and *Metarhizium robertsii*.

The analysis was performed in several consecutive stages.

### 1. Differential gene expression analysis

Differential gene expression was analysed in R using DESeq2.
The analysis was performed separately for haemocytes and fat body. For each tissue, the following comparisons were evaluated:

- *Beauveria bassiana* vs control;
- *Metarhizium robertsii* vs control;
- *B. bassiana* vs *M. robertsii*.

The script used for this analysis is available here: `scripts/01_differential_expression.R`
The input files required to run this script are provided in the `data/` directory:

- `ldec_dge_counts.tsv` — gene-level read count matrix;
- `Design_table.tsv` — experimental design and sample metadata;
- `ldec_gff_annotated.tsv` — gene annotation table.

### 2. Functional annotation

The functional annotation of *L. decemlineata* genes was expanded by integrating information from: NCBI RefSeq annotation, InterProScan, ProteInfer.

InterProScan and ProteInfer predictions were generated from the predicted protein products of the Ldec_2.0 genome assembly (GCF_000500325.1).

Gene identifiers and Gene Ontology terms obtained from the three annotation sources were subsequently merged to generate an extended gene-to-GO annotation dataset.

The input files used at this stage include:

- `GCF_000500325.1_Ldec_2.0_genomic.gff` — NCBI RefSeq genome annotation for the Ldec_2.0 assembly;
- `GCF_000500325.1_Ldec_2.0_gene_ontology.tab` — NCBI Gene Ontology annotation for *L. decemlineata*;
- `gene2go` — NCBI Gene-to-GO association file;
- `go-basic.obo` — Gene Ontology structure file;
- `goslim_generic.obo` — generic GO slim ontology;
- `interpro.csv` — functional annotation generated using InterProScan;
- `protinfer.csv` — functional predictions generated using ProteInfer.

Reference genome annotation and protein sequences for the Ldec_2.0 assembly were obtained from NCBI RefSeq (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/500/325/GCF_000500325.1_Ldec_2.0/)
Gene Ontology resources were obtained from the Gene Ontology Consortium and NCBI Gene (https://geneontology.org/docs/download-ontology/).
