#!/usr/bin/env Rscript

required_packages <- c(
  "iNEXT", "lme4", "lmerTest", "MuMIn", "randomForest", "vegan",
  "ggplot2", "patchwork", "pdp", "corrplot"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) == 0L) {
  message("All required packages are already installed.")
} else {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
