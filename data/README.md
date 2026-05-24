## Data

The `data/` directory contains released processed inputs, compact simulation tables, and links to public source datasets used in the manuscript. Large processed files are managed with Git LFS when needed.

The current release includes processed input and annotation files for the local heme exposure example, together with compact simulation files for the synthetic benchmarks. For other real data analyses, links to the public source datasets are provided, and additional processed examples will be added in future updates.

Please clone the repository with Git LFS enabled before running notebooks that require large processed files:

```bash
git clone https://github.com/liuyk669/AReST2.git
cd AReST2
git lfs pull
