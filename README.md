# Bacterial diversity in European mountain forest habitats

[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/Code%20license-MIT-green.svg)](LICENSE)
[![NCBI BioProject](https://img.shields.io/badge/NCBI-PRJNA1336446-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446)

Reproducible R code and processed analysis inputs for the article **“Bacterial Diversity Patterns and Their Environmental Predictors in European Mountain Forest Habitats.”**

The study compares bacterial alpha diversity in paired forest litter and soil samples collected across 17 European mountain areas. It evaluates habitat differences, latitudinal and elevational patterns, phylum-level responses, environmental correlates, and habitat-specific Random Forest predictors.

## Repository contents

```text
.
├── reproduce_analysis.R
├── install_packages.R
├── data/
│   ├── README.md                        data dictionary and derived variables
│   ├── phyloseq_16S_OTU.rds
│   ├── canopy_vegetation.csv
│   ├── climate_CHELSA_by_plot.csv
│   └── elevation_reference_plots.csv
├── results/                             output of the verified run
│   ├── RUN_SUMMARY.txt
│   ├── sessionInfo.txt
│   ├── package_versions.csv
│   ├── input_checksums.csv
│   ├── figures/
│   └── tables/
├── CITATION.cff
└── LICENSE
```

- `reproduce_analysis.R` runs the complete statistical workflow and creates the main and supplementary figures, result tables, fitted model objects, a run summary, and R session information.
- `data/phyloseq_16S_OTU.rds` contains the processed bacterial OTU table, taxonomy, and sample-level environmental metadata.
- `data/canopy_vegetation.csv` contains the canopy species matrix used for the vegetation PCA.
- `data/elevation_reference_plots.csv` contains the G1 reference plots that define relative elevation.
- `data/climate_CHELSA_by_plot.csv` lists the CHELSA climate values of each analysed plot as a plain table.
- `data/README.md` documents all inputs, including how relative elevation and canopy PCA scores are derived and the source of the CHELSA climate variables.
- Raw sequences are available through [NCBI BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446).

## Analyses reproduced

The script includes:

- derivation of relative elevation from the G1 reference plots, with a table of the divisor of each location;
- rarefied Shannon diversity comparison between litter and soil, repeated with coverage-standardized diversity at 0.99 coverage;
- a coverage scan from 0.90 to 0.99 for the whole community and every phylum–habitat combination, with retained sample counts and the extent of extrapolation;
- Gaussian linear mixed-effects models of latitude and relative elevation (AIC comparison of eight candidate models), Benjamini–Hochberg (BH) correction across the four primary habitat × gradient tests, and sensitivity analyses using raw elevation and 0.90 coverage;
- climate-only, soil-only, and combined residual analyses with relative and raw elevation, reporting estimates, standard errors, p-values, and singular fits;
- phylum-level models for all ten major phyla in both habitats, BH-corrected tests of each gradient and of quadratic (curvature) terms, and residual analyses;
- the canopy vegetation PCA including all tree species;
- Random Forest models with rarefied Shannon diversity as the response, OOB performance, variable importance, and partial dependence; and
- Figs. 1–5 and Figs. S1–S4.

## Requirements

The analysis was verified with R 4.5.2 and the package versions in `results/package_versions.csv` (full details in `results/sessionInfo.txt`). The required packages are `iNEXT` (≥ 3.0.0), `lme4`, `lmerTest`, `randomForest`, `vegan`, `ggplot2`, `patchwork`, `pdp`, and `corrplot`; their current CRAN releases require R ≥ 4.1. `MuMIn` is no longer needed. Install the packages with:

```bash
Rscript install_packages.R
```

The installer also reports any package whose installed version differs from the verified run.

## Run the analysis

```bash
git clone https://github.com/branysot/european-mountain-bacterial-diversity.git
cd european-mountain-bacterial-diversity
Rscript reproduce_analysis.R
```

A complete run takes about two minutes. Outputs are written to `results/`; `results/models/fitted_models.rds` is not committed because of its size.

For a faster test run that omits the partial-dependence panels:

```bash
Rscript reproduce_analysis.R --skip-pdp
```

To select a different output directory:

```bash
Rscript reproduce_analysis.R --output=/path/to/output
```

## Changes in the revised workflow

- The canopy PCA now includes *Betula pendula*, which an indexing error had omitted, and its scores are used directly by the Random Forest models.
- Relative elevation is derived in the script from `data/elevation_reference_plots.csv` instead of being read as a precomputed column.
- Coverage estimates are retained when the achieved coverage equals the target within 1e-4. The earlier tolerance of 1e-8 could drop a sample whose interpolated coverage was 0.98999997, depending on platform.
- All ten phyla are analyzed in both habitats under the same eligibility rule.
- Polynomial terms are tested jointly, and quadratic terms are tested and BH-corrected separately.
- Phylum trend lines in Figs. 2–4 show fixed-effect predictions of the fitted mixed models.

## Methodological scope

This repository starts from the processed OTU, taxonomy, and metadata object. Raw 16S sequence processing (quality filtering, chimera removal, OTU clustering, and taxonomic assignment) was performed in SEED as described in the article and is outside the scope of this R script. Climate values were extracted from CHELSA before this workflow; see `data/README.md`.

The analysis uses fixed random seeds. The script checks the processed dataset dimensions and sequencing depth before fitting models and verifies the derived relative elevation against the stored values.

## Data availability

Raw sequence reads are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446). The compact processed inputs needed to reproduce the statistical analyses are included here.

## Citation

Please cite the associated article when using this code or the processed analysis inputs. Citation metadata are provided in [`CITATION.cff`](CITATION.cff); the journal citation and DOI can be added once assigned.

## License

The analysis code is released under the [MIT License](LICENSE). The included processed data are provided for reproducibility of the associated article; please cite the article and the NCBI BioProject when reusing them.
