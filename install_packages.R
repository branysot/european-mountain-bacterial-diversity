#!/usr/bin/env Rscript

# Installs the CRAN packages required by analysis.R and reports any difference
# from the versions used for the published results
# (results/package_versions.csv; R 4.5.2). The packages require R >= 4.1.

required_packages <- c(
  "iNEXT", "lme4", "lmerTest", "randomForest", "vegan",
  "ggplot2", "patchwork", "pdp", "corrplot"
)

if (getRversion() < "4.1.0") {
  stop("R >= 4.1.0 is required; the published results were produced with R 4.5.2.")
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
  published <- read.csv(versions_file)
  published <- published[published$package %in% required_packages, ]
  installed <- vapply(published$package, function(x) {
    tryCatch(as.character(utils::packageVersion(x)), error = function(e) NA_character_)
  }, character(1))
  different <- is.na(installed) | installed != published$version
  if (any(different)) {
    message("Installed versions that differ from those used for the published results:")
    print(data.frame(package = published$package[different],
                     published = published$version[different],
                     installed = installed[different]), row.names = FALSE)
  } else {
    message("Installed versions match those used for the published results.")
  }
}
