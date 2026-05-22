# inoa-sertoli-developmental-projection

Code and source data for developmental-reference projection analysis identifying an immature-like, stress-biased Sertoli-cell state in idiopathic non-obstructive azoospermia.

## Overview

This repository contains analysis scripts and reproducibility materials for a single-cell transcriptomic reanalysis of public human testicular datasets. The study constructs normal postnatal developmental references across major human testicular cell types and projects disease-derived cells from idiopathic non-obstructive azoospermia onto these reference axes.

The main biological focus is the Sertoli-cell maturation continuum. The analysis evaluates whether Sertoli cells from idiopathic non-obstructive azoospermia occupy earlier positions along a normal postnatal Sertoli-cell developmental reference and whether these projected states are associated with stress-response, senescence-associated, support-related and regulatory-activity programs.

This repository is intended to support peer review and reproducibility of the manuscript.

## Manuscript

Developmental-reference projection reveals an immature-like, stress-biased Sertoli-cell state in idiopathic non-obstructive azoospermia

## Public datasets

The analysis reuses public human testicular single-cell RNA sequencing datasets from the Gene Expression Omnibus:

- GSE149512
- GSE182786

No new raw sequencing data were generated in this study.

## Repository structure

```text
.
├── README.md
├── code/
│   ├── 01_preprocessing_and_annotation.R
│   ├── 02_scenic_reference_regulons.R
│   ├── 03_developmental_reference_construction.R
│   ├── 04_disease_cell_projection.R
│   ├── 05_sertoli_molecular_characterization.R
│   └── 06_integrated_sertoli_evidence.R
├── data_access/
│   └── public_dataset_sources.md
├── metadata/
│   ├── sample_metadata.csv
│   ├── cell_type_markers.csv
│   └── gene_sets.csv
├── results/
│   ├── developmental_reference_summaries/
│   ├── disease_projection_summaries/
│   ├── sertoli_gene_set_scores/
│   ├── differential_expression/
│   └── regulon_activity/
├── figure_source_data/
│   ├── Figure1/
│   ├── Figure2/
│   └── Figure3/
└── reproducibility/
    ├── session_info_R.txt
    ├── session_info_python.txt
    └── package_versions.csv

