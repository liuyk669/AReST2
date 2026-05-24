## Data

The `data/` directory contains released processed inputs and synthetic benchmark inputs used in the manuscript. Large processed files are managed with Git LFS when needed.

The current release includes processed input and annotation files for the local heme exposure example, together with two synthetic benchmark inputs:

```text
c_shaped_region/
two_separated_regions/
```

Each simulation folder contains:

```text
counts.mtx.gz
genes.tsv
barcodes.tsv
meta.csv
dataset_summary.json
```

These are input files only; precomputed scores, detected regions, metric tables, and figures are not included for the simulation benchmarks.

For other real data analyses, links to the public source datasets are provided, and additional processed examples will be added in future updates.

Please clone the repository with Git LFS enabled before running notebooks that require large processed files:

```bash
git clone https://github.com/liuyk669/AReST2.git
cd AReST2
git lfs pull
