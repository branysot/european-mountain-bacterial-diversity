# Processed analysis inputs

This directory contains the compact inputs needed by `reproduce_analysis.R`.

## `phyloseq_16S_OTU.rds`

Processed bacterial community object containing:

- an OTU count matrix;
- SILVA-based taxonomic assignments; and
- sample metadata, including habitat, location, elevation, climate, soil properties, and vegetation PCA scores.

The object contains 318 samples before the filtering and pairing steps applied by the analysis script. It is read directly from its stored components, so the `phyloseq` package is not required.

SHA-256: `e5955364fe5af8bad99437c8a5b6eac9a5b47d2f2540bc78540a7a14ca14bdb2`

## `canopy_vegetation.csv`

Canopy species matrix used to reconstruct the vegetation PCA and species loadings. It was extracted from the `Vegetation_canopy` worksheet of the study metadata workbook; unrelated laboratory and sequencing worksheets are not included.

SHA-256: `deace20fb5d95586899727989be256b68bffa25469258651259ef1993af3ebcd`

Raw sequence reads are available from [NCBI BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446).
