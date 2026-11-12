# Bacterial diversity in European mountain forest habitats

[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/Code%20license-MIT-green.svg)](LICENSE)
[![NCBI BioProject](https://img.shields.io/badge/NCBI-PRJNA1336446-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446)

R code and processed data for the article **“Bacterial Diversity Patterns and Their Environmental Predictors in European Mountain Forest Habitats.”**

The study compares bacterial alpha diversity in paired forest litter and soil samples from 17 European mountain areas. It tests the difference between the two habitats, the patterns along latitude and elevation, the responses of individual phyla, and the environmental predictors of diversity in each habitat.

## Repository contents

```text
.
├── analysis.R             complete statistical analysis and figures
├── install_packages.R     installs the required R packages
├── data/                  processed inputs (described in data/README.md)
│   ├── phyloseq_16S_OTU.rds
│   ├── canopy_vegetation.csv
│   ├── elevation_reference_plots.csv
│   └── climate_CHELSA_by_plot.csv
├── results/               output of analysis.R
│   ├── tables/
│   ├── figures/
│   ├── RUN_SUMMARY.txt
│   ├── sessionInfo.txt
│   ├── package_versions.csv
│   └── input_checksums.csv
├── CITATION.cff
└── LICENSE
```

Raw sequences are available from [NCBI BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446).

## What the script does

1. Calculates relative elevation from the tree-line reference plots (plot elevation divided by the tree-line elevation of the same location).
2. Compares Shannon diversity between litter and soil using rarefied data, and repeats the comparison with coverage-standardized diversity.
3. Estimates coverage-standardized diversity (Hill number of order q = 1, iNEXT) at sample coverages from 0.90 to 0.99.
4. Fits linear mixed-effects models of diversity along latitude and relative elevation, with Benjamini–Hochberg correction, and repeats them with raw elevation and at 0.90 coverage.
5. Tests whether the geographic patterns remain after accounting for climate, soil properties, or both (residual analyses).
6. Repeats the geographic and residual analyses for the ten most abundant phyla in both habitats.
7. Runs a PCA of the canopy tree species composition.
8. Fits Random Forest models of Shannon diversity for each habitat and reports out-of-bag performance, variable importance, and partial dependence.
9. Creates Figs. 1–5 and Figs. S1–S4 of the article.

## Requirements

R ≥ 4.1 and the packages `iNEXT` (≥ 3.0.0), `lme4`, `lmerTest`, `randomForest`, `vegan`, `ggplot2`, `patchwork`, `pdp`, and `corrplot`. Install them with:

```bash
Rscript install_packages.R
```

The results in `results/` were produced with R 4.5.2 and the package versions listed in `results/package_versions.csv` (full details in `results/sessionInfo.txt`).

## Run the analysis

```bash
git clone https://github.com/branysot/european-mountain-bacterial-diversity.git
cd european-mountain-bacterial-diversity
Rscript analysis.R
```

A complete run takes about two minutes and writes its output to `results/`. The analysis uses fixed random seeds, so repeated runs give the same results. Fitted model objects are saved to `results/models/`, which is not tracked in the repository because of its size.

Options:

```bash
Rscript analysis.R --skip-pdp                  # skip the partial-dependence figures (faster)
Rscript analysis.R --output=/path/to/output    # write the output to another directory
```

## Scope

The analysis starts from the processed OTU table, taxonomy, and sample metadata. Raw 16S sequence processing (quality filtering, chimera removal, OTU clustering, and taxonomic assignment) was performed in SEED as described in the article. Climate values were extracted from the CHELSA database beforehand; see `data/README.md`.

## Citation

Please cite the associated article when using this code or the processed data. Citation metadata are provided in [`CITATION.cff`](CITATION.cff); the journal citation and DOI will be added once assigned.

## License

The code is released under the [MIT License](LICENSE). The processed data are provided for reproducibility of the associated article; please cite the article and the NCBI BioProject when reusing them.
