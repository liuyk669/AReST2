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

- `data/`, which contains processed data used in the manuscript, including real spatial transcriptomics data and simulation related files .
- `notebooks/`, which contains reproduction notebooks and helper code for running AReST analyses.
- `results/`, which contains reproduced outputs, including anomaly scores, detected regions, metric summaries, and representative figures.

Additional datasets, benchmark scripts, and reproduction modules will be added in future updates.
