# Processed data

Inputs of `analysis.R`.

## `phyloseq_16S_OTU.rds`

Processed bacterial community data:

- OTU count matrix;
- SILVA-based taxonomy; and
- sample metadata: habitat, plot, location, coordinates, elevation, climate, and soil properties.

The object holds 318 samples. The script excludes three litter samples with fewer than 5,000 sequences together with the soil samples of the same plots, which leaves 312 samples (156 litter and 156 soil) from 156 plots. The object is read from its stored components, so the `phyloseq` package is not required. The sample metadata contain further columns from the study database that are not used in this analysis.

## `canopy_vegetation.csv`

Canopy species matrix: percentage canopy cover of each tree species within approximately 50 m of each plot, from the botanical survey. The script uses it for the canopy PCA.

## `elevation_reference_plots.csv`

Elevation and coordinates of the 50 plots of the uppermost elevation level (G1) of each gradient. These plots were located at the local tree line. They had no forest litter and are therefore not part of the bacterial dataset. The highest G1 plot elevation of each location is the divisor used to calculate relative elevation.

## `climate_CHELSA_by_plot.csv`

The 19 CHELSA bioclimatic variables of each of the 156 analysed plots, with plot, location, coordinates (decimal degrees), and elevation (m a.s.l.). The same values are stored in the sample metadata of `phyloseq_16S_OTU.rds`, from which the script reads them; this table is provided for easy inspection.

## Variables calculated by the script

| Variable | Section of `analysis.R` | Definition |
|---|---|---|
| Relative elevation (`Elevation_relative`) | 3 | Plot elevation (m a.s.l.) divided by the highest elevation of the G1 reference plots of the same location. 0 = sea level; 1 = local tree line. The divisors are listed in `results/tables/02_relative_elevation_divisors.csv`. |
| Canopy PCA scores (`canopy_PC1`–`canopy_PC8`) | 5 | PCA (`vegan::rda`, `scale = TRUE`) of the canopy species matrix. The first eight axes explain 51.4% of the variance. |
| Coverage-standardized diversity (`q1`) | 6 | Hill number of order q = 1 at 0.99 sample coverage (`iNEXT::estimateD`). |

## Climate variables

The 19 bioclimatic variables (columns `Annual_Mean_Temp` to `Precipitation_Coldest_Quarter`; `MAT` and `MAP` duplicate `Annual_Mean_Temp` and `Annual_Precipitation`) are CHELSA climatologies (Karger et al. 2017, *Scientific Data* 4:170122), extracted at the plot coordinates as in Barbi et al. (2025, *New Phytologist*). The values are supplied as extracted; the extraction is not part of `analysis.R`. Plots of the same elevation level usually fall within the same climate grid cell and therefore share climate values.

Units: temperatures in °C; temperature seasonality as the standard deviation of monthly mean temperatures (°C); isothermality in %; precipitation in mm; precipitation seasonality as the coefficient of variation (%).

## Soil variables

`pH` (in water, 1:10), `c` (organic C), `n` (total N), `p_MEL` (Mehlich 3 available P), and `Ca` were measured on the composite soil sample of each plot, as described in the article.
