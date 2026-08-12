# Bacterial diversity in European mountain forest habitats

[![R](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?logo=r)](https://www.r-project.org/)
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
│   ├── phyloseq_16S_OTU.rds
│   └── canopy_vegetation.csv
├── CITATION.cff
└── LICENSE
```

- `reproduce_analysis.R` runs the complete statistical workflow and creates the main and supplementary figures, result tables, fitted model objects, a run summary, and R session information.
- `data/phyloseq_16S_OTU.rds` contains the processed bacterial OTU table, taxonomy, and sample-level environmental metadata used by the analysis.
- `data/canopy_vegetation.csv` contains the canopy species matrix used for the vegetation PCA.
- Raw sequences are not stored in this repository. They are available through [NCBI BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446).

## Analyses reproduced

The script includes:

- rarefied Shannon diversity comparison between litter and soil;
- coverage-standardized Hill diversity of order `q = 1` using iNEXT;
- mixed-effects models of latitude and relative elevation;
- a sensitivity analysis using raw elevation;
- climate-, soil-, and combined-environment residual analyses;
- phylum-level geographic and residual analyses with multiple-testing correction;
- predictor collinearity screening;
- tuned Random Forest models, OOB R² and RMSE, variable importance, and partial-dependence plots;
- canopy vegetation PCA and species-loading plots; and
- reproduction of Figs. 1–5 and Figs. S1–S4.

## Requirements

The analysis was verified with R 4.5.2 and the package versions recorded in the generated `results/sessionInfo.txt`. Install the required packages with:

```bash
Rscript install_packages.R
```

The required R packages are `iNEXT`, `lme4`, `lmerTest`, `MuMIn`, `randomForest`, `vegan`, `ggplot2`, `patchwork`, `pdp`, and `corrplot`.

## Run the analysis

Clone the repository and run the script from any directory:

```bash
git clone https://github.com/branysot/european-mountain-bacterial-diversity.git
cd european-mountain-bacterial-diversity
Rscript reproduce_analysis.R
```

By default, outputs are written to `results/`:

```text
results/
├── figures/
├── tables/
├── models/fitted_models.rds
├── RUN_SUMMARY.txt
└── sessionInfo.txt
```

For a faster test run that omits the 20 partial-dependence panels:

```bash
Rscript reproduce_analysis.R --skip-pdp
```

To select a different output directory:

```bash
Rscript reproduce_analysis.R --output=/path/to/output
```

The default run reproduces the canopy PCA input used for the article. The optional `--corrected-pca` argument includes *Betula pendula* as a transparent sensitivity analysis.

## Methodological scope

This repository starts from the processed OTU, taxonomy, and metadata object. Raw 16S sequence processing—including quality filtering, chimera removal, OTU clustering, and taxonomic assignment—was performed in SEED as described in the article and is outside the scope of this R script.

The analysis uses a fixed random seed (`set.seed(123)`). The script also checks the processed dataset dimensions and sequencing depth before fitting models, helping detect accidental use of a different input object.

## Data availability

Raw sequence reads are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446). The compact processed inputs needed to reproduce the statistical analyses are included here.

## Citation

Please cite the associated article when using this code or the processed analysis inputs. Citation metadata are provided in [`CITATION.cff`](CITATION.cff); the journal citation and DOI can be added once assigned.

## License

The analysis code is released under the [MIT License](LICENSE). The included processed data are provided for reproducibility of the associated article; please cite the article and the NCBI BioProject when reusing them.
