## Data sources and coordinate alignment

The `data/` directory contains processed data, compact simulation files, and links to public datasets used in the manuscript. Large processed files are managed with Git LFS when needed.

The real spatial transcriptomics datasets used in the manuscript are available from the following public sources:

- Local heme exposure mouse brain dataset: [GSE182127](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182127)
- Spinal cord injury and neurorehabilitation dataset: [GSE184369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184369)
- Mouse cortex stab wound injury dataset: [GSE226208](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226208)

For analyses involving multiple tissue sections or biological replicates, spatial coordinates need to be interpreted in a comparable coordinate system. In this study, tissue sections were aligned to a common spatial coordinate frame before applying AReST. We used STalign, a Python tool for spatial transcriptomics alignment based on diffeomorphic metric mapping:

- STalign: https://github.com/JEFworks-Lab/STalign

Additional processed inputs, benchmark scripts, and reproduction modules will be added in future updates.
