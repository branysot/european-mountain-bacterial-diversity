#!/usr/bin/env Rscript

# Installs the CRAN packages required by reproduce_analysis.R and reports any
# difference from the versions used for the verified run
# (results/package_versions.csv; R 4.5.2).
#
# MuMIn is no longer required: marginal and conditional R2 are computed
# directly in reproduce_analysis.R. The remaining packages declare R >= 4.1 as
# their minimum version.

required_packages <- c(
  "iNEXT", "lme4", "lmerTest", "randomForest", "vegan",
  "ggplot2", "patchwork", "pdp", "corrplot"
)

if (getRversion() < "4.1.0") {
  stop("R >= 4.1.0 is required; the analysis was verified with R 4.5.2.")
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) == 0L) {
  message("All required packages are already installed.")
} else {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
versions_file <- file.path(script_dir, "results", "package_versions.csv")
if (file.exists(versions_file)) {
  verified <- read.csv(versions_file)
  verified <- verified[verified$package %in% required_packages, ]
  installed <- vapply(verified$package, function(x) {
    tryCatch(as.character(utils::packageVersion(x)), error = function(e) NA_character_)
  }, character(1))
  different <- is.na(installed) | installed != verified$version
  if (any(different)) {
    message("Installed versions that differ from the verified run:")
    print(data.frame(package = verified$package[different],
                     verified = verified$version[different],
                     installed = installed[different]), row.names = FALSE)
  } else {
    message("Installed versions match the verified run.")
  }
}
