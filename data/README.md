# Processed analysis inputs

This directory contains the compact inputs needed by `reproduce_analysis.R`.

## `phyloseq_16S_OTU.rds`

Processed bacterial community object containing:

- an OTU count matrix;
- SILVA-based taxonomic assignments; and
- sample metadata, including habitat, location, elevation, climate, and soil properties.

The object contains 318 samples before the filtering and pairing steps applied by the analysis script. It is read directly from its stored components, so the `phyloseq` package is not required.

## `canopy_vegetation.csv`

Canopy species matrix (percentage canopy cover of each tree species within approximately 50 m of each plot) from the botanical survey. It was extracted from the `Vegetation_canopy` worksheet of the study metadata workbook. The script computes the canopy PCA from this file and uses the resulting scores in the Random Forest models. Rows of the G1 reference plots (no forest canopy) are excluded by the script.

## `elevation_reference_plots.csv`

Elevation and coordinates of the 50 G1 reference plots, the uppermost elevation level of each gradient. These plots were located at the local tree line and classified as above the tree line in the plot metadata; they had no forest litter and are therefore not part of the bacterial dataset. The highest G1 plot elevation of each location is the divisor used to compute relative elevation (see below and `results/tables/02_relative_elevation_divisors.csv`).

## `climate_CHELSA_by_plot.csv`

The 19 CHELSA bioclimatic variables of each of the 156 analysed plots, with plot, location, coordinates (decimal degrees), and elevation (m a.s.l.). The values are the same as those stored in the sample metadata of `phyloseq_16S_OTU.rds` (the litter and soil samples of a plot share them) and are provided as a plain table for easy inspection; the script reads them from the `.rds` file. Variables and units are described under "Climate variables" below.

## Derived variables

| Variable | Where it is derived | Definition |
|---|---|---|
| Relative elevation (`Elevation_relative`) | `reproduce_analysis.R`, section 3 | Plot elevation (m a.s.l.) divided by the highest elevation of the G1 reference plots of the same location. 0 = sea level; 1 = highest G1 reference plot. The script verifies that the result equals the `Elevation3` column stored in the metadata. |
| Canopy PCA scores (`canopy_PC1`–`canopy_PC8`) | `reproduce_analysis.R`, section 5 | PCA (`vegan::rda`, `scale = TRUE`) of the canopy species matrix, all species included. The first eight axes explain 51.4% of the variance. |
| Coverage-standardized diversity (`q1`) | `reproduce_analysis.R`, section 6 | Hill number of order q = 1 at 0.99 sample coverage (`iNEXT::estimateD`). |

The metadata also contain the canopy scores used in the first submission (`PC1.x`–`PC8.x`). They were computed with *Betula pendula* inadvertently omitted and are used only to show the effect of the correction (results tables 07 and 22).

## Climate variables

The 19 bioclimatic variables (columns `Annual_Mean_Temp` to `Precipitation_Coldest_Quarter`; `MAT` and `MAP` duplicate `Annual_Mean_Temp` and `Annual_Precipitation`) are CHELSA climatologies (Karger et al. 2017, *Scientific Data* 4:170122), extracted at the plot coordinates, as in Barbi et al. (2025, *New Phytologist*). Units: temperatures in °C, temperature seasonality as the standard deviation of monthly mean temperatures (°C), isothermality in %, precipitation in mm, and precipitation seasonality as the coefficient of variation (%). Because plots of the same elevation level usually fall within the same 30 arc-second grid cell, they share climate values. The values are supplied as extracted (see also `climate_CHELSA_by_plot.csv`); raster extraction is not part of this workflow.

## Soil variables

`pH` (water, 1:10), `c` (organic C), `n` (total N), `p_MEL` (Mehlich 3 available P) and `Ca` were measured on the composite soil sample of each plot, as described in the article.

Raw sequence reads are available from [NCBI BioProject PRJNA1336446](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1336446).
