# AReST

AReST is a spatially weighted neighborhood contrast framework for detecting anomalous regions in comparative spatial transcriptomics data.

Comparative spatial transcriptomics is often used to ask not only which genes change between biological conditions, but also where condition related molecular changes appear within tissue. AReST addresses this problem by estimating spot level anomaly evidence from transcriptional neighborhoods, incorporating spatial proximity, correcting for the global condition prior, and converting positive anomaly evidence into spatially coherent anomalous regions.

## Overview

AReST takes as input:

- a gene expression matrix,
- spatial coordinates,
- condition labels.

The workflow consists of four main steps:

1. Construct transcriptional neighborhoods in expression space.
2. Compute a spatially weighted local condition enrichment score.
3. Derive an upper threshold using label permutation.
4. Apply DBSCAN to significant anomalous spots to obtain anomalous regions.

The detected regions can then be used for downstream biological interpretation, including marker analysis, region level enrichment, and comparison with anatomical or literature supported reference regions.

## Repository organization

This repository is organized into three main parts:

- `data/`, which contains processed data used in the manuscript, including real spatial transcriptomics data and simulation related files.
- `notebooks/`, which contains reproduction notebooks and helper code for running AReST analyses.
- `results/`, which contains reproduced outputs, including anomaly scores, detected regions, metric summaries, and representative figures.

The current release includes processed input and annotation files for the local heme exposure example, together with two synthetic benchmark inputs.

## Simulation Inputs

We considered two representative settings: a single non-convex C-shaped anomalous region and a two-region setting with spatially separated anomalous regions of different perturbation strengths. The corresponding inputs are included as separate folders:

```text
data/c_shaped_region/
data/two_separated_regions/
```

Each simulation folder contains only the files needed to rerun methods, without precomputed scores or output figures:

```text
counts.mtx.gz
genes.tsv
barcodes.tsv
meta.csv
dataset_summary.json
```

`counts.mtx.gz` is a genes-by-spots Matrix Market count matrix. `meta.csv` contains spot barcodes, spatial coordinates, condition labels, replicate labels, and the ground-truth anomaly label.

## Git LFS

The processed `.h5ad` file is stored with Git LFS. Please clone the repository with Git LFS enabled before running notebooks that require large processed files:

```bash
git clone https://github.com/liuyk669/AReST2.git
cd AReST2
git lfs pull
```

The compressed simulation matrices are stored directly in git and do not require Git LFS.
