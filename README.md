# AReST2

Minimal reproducibility materials for the GSE182127 AReST example.

The notebook recomputes AReST anomaly scores from the released `.h5ad` file,
applies permutation max-statistic thresholding, extracts the anomalous region
with DBSCAN, and reports the score-level and region-level metrics used in the
paper example.

## Repository layout

```text
data/
  GSE182127_all_samples_aligned_to_GSM5519060_sham_combined.h5ad
  author_leiden_clusters_heme1000.csv
notebooks/
  GSE182127_AReST_reproduction.ipynb
  gse182127.py
results/
  gse182127_reproduction/
```

The `.h5ad` file is stored with Git LFS. If it is not downloaded
automatically, install Git LFS and run:

```bash
git lfs pull
```

## Run

```bash
pip install -r requirements.txt
jupyter notebook notebooks/GSE182127_AReST_reproduction.ipynb
```

The notebook uses repository-relative paths, so it can be run from a cloned
checkout without editing absolute server paths.

## Main settings

```text
K_expr = 80
Gaussian bandwidth quantile = 0.60
permutations = 1000
alpha = 0.05
DBSCAN eps = median sixth-nearest-neighbor distance among significant spots
DBSCAN min_samples = 7
minimum retained cluster size = 20
```

The reference region is defined from the author-provided Leiden clusters 10
and 12 in `author_leiden_clusters_heme1000.csv`.

## Expected output

The executed notebook should reproduce the included result files under
`results/gse182127_reproduction/`. The key values are:

```text
significant spots = 1121
DBSCAN region size = 898
overlap with Leiden 10/12 reference = 670 / 704
precision = 0.7461
recall = 0.9517
Jaccard index = 0.7189
AUROC = 0.9722
AUPRC = 0.9229
Moran's I = 0.9368
```
