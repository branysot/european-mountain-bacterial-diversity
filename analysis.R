#!/usr/bin/env Rscript

# Analysis code for:
# "Bacterial Diversity Patterns and Their Environmental Predictors in
#  European Mountain Forest Habitats"
#
# Starting point
# --------------
# This script runs the statistical analyses and creates the figures from the
# processed OTU table and associated environmental metadata in `data/`.
# Raw sequence processing was performed in SEED as described in the article
# and is outside the scope of this script.
#
# Run from any directory:
#   Rscript analysis.R
#
# Optional arguments:
#   --output=/path/to/output_directory
#   --skip-pdp       Skip the partial-dependence panels (faster test run)
#
# Required R packages (see install_packages.R and package_versions.csv):
#   iNEXT (>= 3.0.0), lme4, lmerTest, randomForest, vegan, ggplot2,
#   patchwork, pdp, corrplot
#
# Derived variables computed in this script (see data/README.md):
#   - relative elevation (plot elevation / highest G1 reference plot elevation
#     of the same location; data/elevation_reference_plots.csv)
#   - canopy PCA scores (data/canopy_vegetation.csv)
# Climate variables are the CHELSA values stored in the sample metadata.

options(stringsAsFactors = FALSE)
options(contrasts = c("contr.treatment", "contr.poly"))
set.seed(123)

# -----------------------------------------------------------------------------
# 0. Paths, arguments, packages
# -----------------------------------------------------------------------------

command_args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
if (length(script_arg) == 0L) {
  stop("Run this file with Rscript so its data paths can be resolved.")
}

script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
script_dir <- dirname(script_file)

output_arg <- grep("^--output=", command_args, value = TRUE)
output_dir <- if (length(output_arg) == 1L) {
  normalizePath(sub("^--output=", "", output_arg), mustWork = FALSE)
} else {
  file.path(script_dir, "results")
}

run_pdp <- !"--skip-pdp" %in% command_args

data_dir <- file.path(script_dir, "data")
phyloseq_file <- file.path(data_dir, "phyloseq_16S_OTU.rds")
canopy_file <- file.path(data_dir, "canopy_vegetation.csv")
elevation_reference_file <- file.path(data_dir, "elevation_reference_plots.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "models")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "iNEXT", "lme4", "lmerTest", "randomForest", "vegan",
  "ggplot2", "patchwork", "pdp", "corrplot"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "), "\n",
    "Install them with: Rscript install_packages.R"
  )
}
if (utils::packageVersion("iNEXT") < "3.0.0") {
  stop("iNEXT >= 3.0.0 is required (installed: ",
       as.character(utils::packageVersion("iNEXT")), ").")
}

suppressPackageStartupMessages({
  library(iNEXT)
  library(lme4)
  library(lmerTest)
  library(randomForest)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(pdp)
  library(corrplot)
})

message("Output directory: ", output_dir)

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------

write_table <- function(x, filename, row.names = FALSE) {
  write.csv(x, file.path(table_dir, filename), row.names = row.names,
            na = "NA", fileEncoding = "UTF-8")
}

strip_stored_class <- function(x) {
  # A phyloseq object is an S4 object whose slots are stored as attributes.
  # Removing the stored component class permits extraction even when phyloseq
  # is not installed; this substantially reduces the dependencies of this file.
  attr(x, "class") <- NULL
  x
}

read_phyloseq_components <- function(filename) {
  object <- readRDS(filename)
  otu_object <- attr(object, "otu_table")
  taxa_are_rows <- isTRUE(attr(otu_object, "taxa_are_rows"))
  otu <- strip_stored_class(otu_object)
  tax <- strip_stored_class(attr(object, "tax_table"))
  metadata <- as.data.frame(
    strip_stored_class(attr(object, "sam_data")),
    stringsAsFactors = FALSE
  )

  if (taxa_are_rows) otu <- t(otu)
  storage.mode(otu) <- "numeric"
  tax <- as.matrix(tax)

  if (!identical(rownames(metadata), rownames(otu))) {
    # The stored sample_data component uses sequential row names, while its
    # `Samples` column contains the sequencing identifiers used by otu_table.
    id_candidates <- names(metadata)[vapply(metadata, function(x) {
      length(intersect(as.character(x), rownames(otu))) == nrow(otu)
    }, logical(1))]
    if (length(id_candidates) > 0L) {
      rownames(metadata) <- as.character(metadata[[id_candidates[[1]]]])
    } else if (nrow(metadata) == nrow(otu)) {
      # phyloseq guarantees aligned component order; use it as a final fallback.
      rownames(metadata) <- rownames(otu)
    }
    common_samples <- intersect(rownames(otu), rownames(metadata))
    otu <- otu[common_samples, , drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]
  }
  if (!identical(colnames(otu), rownames(tax))) {
    common_taxa <- intersect(colnames(otu), rownames(tax))
    otu <- otu[, common_taxa, drop = FALSE]
    tax <- tax[common_taxa, , drop = FALSE]
  }

  list(otu = otu, tax = tax, metadata = metadata)
}

shannon_index <- function(count_matrix) {
  proportions <- count_matrix / rowSums(count_matrix)
  -rowSums(ifelse(proportions > 0, proportions * log(proportions), 0))
}

clean_phylum_names <- function(x) {
  x <- sub("^p__", "", x)
  mapping <- c(
    Proteobacteria = "Pseudomonadota",
    Actinobacteriota = "Actinomycetota",
    Chloroflexi = "Chloroflexota",
    Patescibacteria = "Patescibacteriota"
  )
  replace <- !is.na(x) & x %in% names(mapping)
  x[replace] <- unname(mapping[x[replace]])
  x
}

# Marginal and conditional R2 for Gaussian LMMs (Nakagawa & Schielzeth 2013):
# variance of the fixed-effect predictions, the summed random-intercept
# variances, and the residual variance. For the random-intercept models used
# here this equals the output of MuMIn::r.squaredGLMM().
lmm_r2 <- function(model) {
  fixed_prediction <- as.vector(lme4::getME(model, "X") %*% lme4::fixef(model))
  variance_fixed <- stats::var(fixed_prediction)
  components <- as.data.frame(lme4::VarCorr(model))
  variance_random <- sum(components$vcov[components$grp != "Residual"])
  variance_residual <- stats::sigma(model)^2
  total <- variance_fixed + variance_random + variance_residual
  c(marginal = variance_fixed / total,
    conditional = (variance_fixed + variance_random) / total)
}

# Test of one gradient term. For a linear term, the F-test equals the t-test of
# its coefficient. For an orthogonal quadratic polynomial, the F-test jointly
# tests both polynomial coefficients; the quadratic coefficient is also
# returned so that curvature can be evaluated separately.
term_test <- function(model, variable) {
  term_table <- suppressMessages(anova(model, type = 3))
  coefficients <- coef(summary(model))
  p_column <- grep("^Pr\\(", colnames(coefficients), value = TRUE)[[1]]
  polynomial_term <- paste0("poly(", variable, ", 2)")
  result <- list(term = NA_character_, polynomial = NA, estimate = NA_real_,
                 se = NA_real_, p = NA_real_, quadratic_estimate = NA_real_,
                 quadratic_se = NA_real_, quadratic_p = NA_real_)
  if (variable %in% rownames(term_table)) {
    result$term <- variable
    result$polynomial <- FALSE
    result$estimate <- coefficients[variable, "Estimate"]
    result$se <- coefficients[variable, "Std. Error"]
    result$p <- term_table[variable, "Pr(>F)"]
  } else if (polynomial_term %in% rownames(term_table)) {
    quadratic_row <- paste0(polynomial_term, "2")
    result$term <- polynomial_term
    result$polynomial <- TRUE
    result$p <- term_table[polynomial_term, "Pr(>F)"]
    result$quadratic_estimate <- coefficients[quadratic_row, "Estimate"]
    result$quadratic_se <- coefficients[quadratic_row, "Std. Error"]
    result$quadratic_p <- coefficients[quadratic_row, p_column]
  }
  result
}

interaction_test <- function(model) {
  term_table <- suppressMessages(anova(model, type = 3))
  rows <- grepl(":", rownames(term_table), fixed = TRUE)
  if (!any(rows)) return(NA_real_)
  min(term_table[rows, "Pr(>F)"], na.rm = TRUE)
}

candidate_formulas <- function(response, elevation_variable, interaction = FALSE) {
  operator <- if (interaction) "*" else "+"
  random <- "+ (1 | PermGroup) + (1 | Location)"
  list(
    M1 = as.formula(sprintf(
      "%s ~ poly(%s, 2) %s poly(Latitude, 2) %s",
      response, elevation_variable, operator, random
    )),
    M2 = as.formula(sprintf(
      "%s ~ poly(%s, 2) %s Latitude %s",
      response, elevation_variable, operator, random
    )),
    M3 = as.formula(sprintf(
      "%s ~ %s %s poly(Latitude, 2) %s",
      response, elevation_variable, operator, random
    )),
    M4 = as.formula(sprintf(
      "%s ~ %s %s Latitude %s",
      response, elevation_variable, operator, random
    ))
  )
}

# Candidate models are compared by AIC under maximum likelihood; the selected
# model is refitted with REML for inference. PermGroup is the elevation level
# within a location, so the random intercepts represent elevation level nested
# within location.
fit_candidate_lmms <- function(data, response = "sqrt_q1",
                               elevation_variable = "Elevation_relative",
                               interaction = FALSE) {
  formulas <- candidate_formulas(response, elevation_variable, interaction)
  ml_models <- lapply(formulas, function(formula) {
    suppressMessages(suppressWarnings(lmer(formula, data = data, REML = FALSE)))
  })
  aic_values <- vapply(ml_models, AIC, numeric(1))
  best_name <- names(which.min(aic_values))
  best_reml <- suppressMessages(suppressWarnings(lmer(
    formulas[[best_name]], data = data, REML = TRUE
  )))

  list(
    formulas = formulas,
    ml_models = ml_models,
    aic = aic_values,
    best_name = best_name,
    best_model = best_reml,
    elevation = term_test(best_reml, elevation_variable),
    latitude = term_test(best_reml, "Latitude"),
    interaction_p = if (interaction) interaction_test(best_reml) else NA_real_,
    singular = lme4::isSingular(best_reml, tol = 1e-4),
    r2 = lmm_r2(best_reml)
  )
}

candidate_aic_table <- function(fit, habitat, elevation_type, interaction) {
  data.frame(
    habitat = habitat,
    elevation_type = elevation_type,
    interaction = interaction,
    model = names(fit$aic),
    formula = vapply(fit$formulas, function(x) paste(deparse(x), collapse = " "),
                     character(1)),
    AIC = unname(fit$aic),
    selected_within_set = names(fit$aic) == fit$best_name
  )
}

gradient_columns <- function(test, prefix) {
  out <- data.frame(
    test$term, test$polynomial, test$estimate, test$se, test$p,
    test$quadratic_estimate, test$quadratic_se, test$quadratic_p
  )
  names(out) <- paste0(prefix, c("_term", "_polynomial", "_estimate", "_SE",
                                 "_p", "_quadratic_estimate",
                                 "_quadratic_SE", "_quadratic_p"))
  out
}

model_summary_row <- function(fit, habitat, elevation_type, interaction,
                              n_samples) {
  cbind(
    data.frame(
      habitat = habitat,
      elevation_type = elevation_type,
      interaction = interaction,
      n_samples = n_samples,
      selected_model = fit$best_name,
      AIC = min(fit$aic)
    ),
    gradient_columns(fit$elevation, "elevation"),
    gradient_columns(fit$latitude, "latitude"),
    data.frame(
      interaction_p = fit$interaction_p,
      singular_fit = fit$singular,
      R2_marginal = fit$r2[["marginal"]],
      R2_conditional = fit$r2[["conditional"]]
    )
  )
}

# Coverage-based q = 1 Hill numbers. iNEXT computes the sample size m at which
# each sample reaches the target coverage: by interpolation when the observed
# coverage exceeds the target and by extrapolation otherwise. Estimates are
# retained when they are finite and the achieved coverage equals the target
# within `sc_tolerance`; the tolerance accommodates the numerical optimization
# iNEXT uses for interpolation (achieved values such as 0.98999997), which can
# otherwise drop samples in a platform-dependent manner.
coverage_estimates <- function(sample_by_taxon, levels, sc_tolerance = 1e-4) {
  taxon_by_sample <- t(sample_by_taxon)
  taxon_by_sample <- taxon_by_sample[rowSums(taxon_by_sample) > 1, , drop = FALSE]
  taxon_by_sample <- taxon_by_sample[, colSums(taxon_by_sample) > 1, drop = FALSE]
  if (nrow(taxon_by_sample) < 2L || ncol(taxon_by_sample) < 2L) return(NULL)
  estimates <- suppressWarnings(iNEXT::estimateD(
    taxon_by_sample, datatype = "abundance", base = "coverage",
    level = levels, nboot = 0
  ))
  estimates <- estimates[estimates$Order.q == 1, , drop = FALSE]
  # Within each assemblage, iNEXT returns the coverage levels in input order.
  level_index <- ave(seq_len(nrow(estimates)), estimates$Assemblage, FUN = seq_along)
  stopifnot(max(level_index) == length(levels))
  estimates$goalSC <- levels[level_index]
  information <- iNEXT::DataInfo(taxon_by_sample, datatype = "abundance")
  reads <- setNames(information$n, information$Assemblage)
  observed_coverage <- setNames(information$SC, information$Assemblage)
  data.frame(
    sample_id = as.character(estimates$Assemblage),
    target_coverage = estimates$goalSC,
    achieved_coverage = estimates$SC,
    reads = unname(reads[as.character(estimates$Assemblage)]),
    observed_coverage = unname(observed_coverage[as.character(estimates$Assemblage)]),
    size_m = estimates$m,
    method = estimates$Method,
    q1 = estimates$qD,
    retained = is.finite(estimates$qD) & is.finite(estimates$m) &
      abs(estimates$SC - estimates$goalSC) <= sc_tolerance,
    n_otus_in_subset = nrow(taxon_by_sample),
    n_libraries_in_subset = ncol(taxon_by_sample)
  )
}

estimate_coverage_q1 <- function(sample_by_taxon, metadata, coverage = 0.99) {
  estimates <- coverage_estimates(sample_by_taxon, levels = coverage)
  estimates <- estimates[estimates$retained, , drop = FALSE]
  sample_ids <- intersect(estimates$sample_id, rownames(metadata))
  estimates <- estimates[match(sample_ids, estimates$sample_id), , drop = FALSE]
  result <- metadata[sample_ids, , drop = FALSE]
  result$sample_id <- sample_ids
  result$q1 <- estimates$q1
  result$sqrt_q1 <- sqrt(estimates$q1)
  result$reads_in_subset <- estimates$reads
  result$observed_coverage <- estimates$observed_coverage
  result$size_m <- estimates$size_m
  result$extrapolation_factor <- estimates$size_m / estimates$reads
  result
}

format_p <- function(x, prefix = "p") {
  ifelse(is.na(x), paste(prefix, "= NA"),
         ifelse(x < 0.001, paste(prefix, "< 0.001"),
                sprintf("%s = %.3f", prefix, x)))
}

theme_figures <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# Fixed-effect prediction of a selected LMM along one gradient, with the other
# gradient held at its mean, back-transformed from sqrt(q1) to q1. The 95%
# confidence band uses the fixed-effect covariance matrix. These lines replace
# simple regressions of q1 on a single gradient, which ignore the other
# gradient and the random effects and can therefore differ in direction from
# the fitted model.
lmm_trend <- function(model, data, gradient, n_points = 100) {
  fixed_terms <- delete.response(terms(model, fixed.only = TRUE))
  variables <- all.vars(fixed_terms)
  grid <- data.frame(seq(min(data[[gradient]]), max(data[[gradient]]),
                         length.out = n_points))
  names(grid) <- gradient
  for (variable in setdiff(variables, gradient)) grid[[variable]] <- mean(data[[variable]])
  design <- model.matrix(fixed_terms, model.frame(fixed_terms, grid))
  fit <- as.vector(design %*% lme4::fixef(model))
  se <- sqrt(rowSums((design %*% as.matrix(stats::vcov(model))) * design))
  stopifnot(isTRUE(all.equal(fit, unname(predict(model, newdata = grid, re.form = NA)))))
  data.frame(
    x = grid[[gradient]],
    q1 = fit^2,
    lower = pmax(fit - 1.96 * se, 0)^2,
    upper = (fit + 1.96 * se)^2
  )
}

# -----------------------------------------------------------------------------
# 2. Load and filter the final bacterial OTU dataset
# -----------------------------------------------------------------------------

if (!file.exists(phyloseq_file)) stop("Missing input: ", phyloseq_file)
components <- read_phyloseq_components(phyloseq_file)
otu <- components$otu
tax <- components$tax
metadata <- components$metadata

keep_bacterial_taxa <-
  (tax[, "Domain"] == "d__Bacteria" | is.na(tax[, "Domain"])) &
  (tax[, "Phylum"] != "p__Chloroplast" | is.na(tax[, "Phylum"])) &
  (tax[, "Phylum"] != "p__Mitochondria" | is.na(tax[, "Phylum"]))
otu <- otu[, keep_bacterial_taxa, drop = FALSE]
tax <- tax[keep_bacterial_taxa, , drop = FALSE]

# Three litter libraries had <=5,000 reads. Their paired soil samples are also
# removed to maintain the paired sampling design.
excluded_plots <- c("FEL3C", "TREL2A", "DAEL3B")
sample_depth_before_pairing <- rowSums(otu)
keep_samples <- sample_depth_before_pairing > 5000 &
  !(metadata$Plot %in% excluded_plots)
otu <- otu[keep_samples, , drop = FALSE]
metadata <- metadata[keep_samples, , drop = FALSE]

# Remove taxa represented by a single retained read after sample filtering.
keep_non_singletons <- colSums(otu) > 1
otu <- otu[, keep_non_singletons, drop = FALSE]
tax <- tax[keep_non_singletons, , drop = FALSE]

metadata$habitat <- factor(metadata$Type.x, levels = c("litter", "soil"))
metadata$Location <- factor(metadata$Location)
metadata$PermGroup <- factor(sub(".$", "", metadata$Plot))

dataset_summary <- data.frame(
  statistic = c(
    "samples", "litter_samples", "soil_samples", "plots", "taxa",
    "total_reads", "minimum_reads", "maximum_reads", "mean_reads"
  ),
  value = c(
    nrow(otu), sum(metadata$habitat == "litter"),
    sum(metadata$habitat == "soil"), length(unique(metadata$Plot)), ncol(otu),
    sum(otu), min(rowSums(otu)), max(rowSums(otu)), mean(rowSums(otu))
  )
)
write_table(dataset_summary, "01_dataset_summary.csv")

stopifnot(
  nrow(otu) == 312L,
  sum(metadata$habitat == "litter") == 156L,
  sum(metadata$habitat == "soil") == 156L,
  sum(otu) == 3993630,
  min(rowSums(otu)) == 5061
)
message("Dataset validation passed: 312 samples and 3,993,630 retained reads.")

# -----------------------------------------------------------------------------
# 3. Relative elevation
# -----------------------------------------------------------------------------
# The uppermost elevation level (G1) of each gradient was located at the local
# tree line and classified as above the tree line in the plot metadata. These
# plots had no forest litter and are not part of the bacterial dataset; the
# highest G1 plot elevation of each location is the divisor that defines
# relative elevation (1 = highest G1 reference plot of that location).

if (!file.exists(elevation_reference_file)) {
  stop("Missing elevation reference input: ", elevation_reference_file)
}
elevation_reference <- read.csv(elevation_reference_file)
stopifnot(all(elevation_reference$Elevation_level == "G1"))
divisors <- tapply(elevation_reference$Elevation, elevation_reference$Location, max)
missing_locations <- setdiff(levels(metadata$Location), names(divisors))
if (length(missing_locations) > 0L) {
  stop("No G1 reference elevation for: ", paste(missing_locations, collapse = ", "))
}
metadata$elevation_divisor <- as.numeric(divisors[as.character(metadata$Location)])
metadata$Elevation_relative <- metadata$Elevation / metadata$elevation_divisor


plot_metadata <- metadata[!duplicated(metadata$Plot), , drop = FALSE]
divisor_table <- do.call(rbind, lapply(split(plot_metadata, plot_metadata$Location),
                                       function(x) {
  location <- as.character(x$Location[1])
  reference <- elevation_reference[elevation_reference$Location == location, ]
  data.frame(
    location = location,
    mountain_range = x$Mountain_range[1],
    latitude = round(mean(x$Latitude), 2),
    n_G1_reference_plots = nrow(reference),
    G1_reference_elevation_range = paste(range(reference$Elevation), collapse = "-"),
    divisor_m = max(reference$Elevation),
    n_forest_plots_analysed = nrow(x),
    forest_plot_elevation_range = paste(range(x$Elevation), collapse = "-"),
    relative_elevation_range = paste(sprintf("%.3f", range(x$Elevation_relative)),
                                     collapse = "-")
  )
}))
write_table(divisor_table, "02_relative_elevation_divisors.csv")

# -----------------------------------------------------------------------------
# 4. Figure 1: rarefied Shannon diversity by habitat
# -----------------------------------------------------------------------------

set.seed(123)
rarefaction_depth <- 5061L
otu_rarefied <- vegan::rrarefy(otu, sample = rarefaction_depth)
shannon_rarefied <- vegan::diversity(otu_rarefied, index = "shannon")

habitat_data <- data.frame(
  sample_id = rownames(metadata),
  shannon = shannon_rarefied[rownames(metadata)],
  habitat = metadata$habitat,
  plot = metadata$Plot,
  location = metadata$Location
)

habitat_model <- lmer(shannon ~ habitat + (1 | location),
                      data = habitat_data, REML = TRUE)
habitat_anova <- anova(habitat_model)
habitat_p <- habitat_anova[["Pr(>F)"]][1]

habitat_statistics <- do.call(rbind, lapply(split(habitat_data, habitat_data$habitat),
                                            function(x) {
  data.frame(
    habitat = as.character(x$habitat[1]), n = nrow(x),
    mean = mean(x$shannon), sd = sd(x$shannon),
    median = median(x$shannon), minimum = min(x$shannon), maximum = max(x$shannon)
  )
}))
habitat_statistics$LMM_soil_minus_litter <- fixef(habitat_model)[["habitatsoil"]]
habitat_statistics$LMM_p <- habitat_p
write_table(habitat_statistics, "03_habitat_shannon_statistics.csv")

habitat_labels <- data.frame(
  habitat = factor(c("litter", "soil"), levels = levels(habitat_data$habitat)),
  label = c("b", "a"),
  y = vapply(split(habitat_data$shannon, habitat_data$habitat), max, numeric(1)) + 0.25
)

fig1 <- ggplot(habitat_data, aes(habitat, shannon, fill = habitat)) +
  geom_boxplot(width = 0.55, outlier.shape = 16) +
  geom_jitter(width = 0.16, alpha = 0.35, size = 1.4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.16,
               linewidth = 0.6, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
               fill = "white", color = "black") +
  geom_text(data = habitat_labels, aes(x = habitat, y = y, label = label),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c(litter = "#01665e", soil = "#8c510a")) +
  coord_cartesian(ylim = c(NA, max(habitat_labels$y) + 0.08)) +
  labs(x = "Habitat", y = "Bacterial diversity (Shannon)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(figure_dir, "Fig1.pdf"), fig1,
       width = 6, height = 5, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig1.png"), fig1,
       width = 6, height = 5, dpi = 220)

# -----------------------------------------------------------------------------
# 5. Canopy vegetation PCA (used by the Random Forest analysis)
# -----------------------------------------------------------------------------
# PCA of the canopy species matrix (all tree species, scaled to unit variance).
# The first eight axes are used as predictors in the Random Forest models.

if (!file.exists(canopy_file)) stop("Missing canopy input: ", canopy_file)
canopy_raw <- read.csv(canopy_file, check.names = FALSE,
                       stringsAsFactors = FALSE)
non_species_columns <- c("ID_sample", "Altitude", "Plot")
blank_columns <- grep("^\\.\\.\\.", names(canopy_raw), value = TRUE)
species_columns <- setdiff(names(canopy_raw), c(non_species_columns, blank_columns))
canopy_species <- as.data.frame(
  lapply(canopy_raw[, species_columns, drop = FALSE],
         function(x) suppressWarnings(as.numeric(x))),
  check.names = FALSE
)
rownames(canopy_species) <- canopy_raw$ID_sample
# G1 reference plots (identifiers with elevation level 1) have no forest canopy.
canopy_species <- canopy_species[!grepl("1", rownames(canopy_species)), , drop = FALSE]
canopy_species <- canopy_species[, colSums(canopy_species, na.rm = TRUE) > 0,
                                 drop = FALSE]
canopy_matrix <- as.matrix(canopy_species)

canopy_pca <- vegan::rda(canopy_matrix, scale = TRUE)
canopy_scores <- as.data.frame(scores(canopy_pca, display = "sites", choices = 1:8))
canopy_loadings <- as.data.frame(scores(canopy_pca, display = "species", choices = 1:8))
canopy_variance <- summary(canopy_pca)$cont$importance[
  "Proportion Explained", 1:8
]

canopy_variance_table <- data.frame(
  axis = paste0("PC", 1:8),
  proportion_explained = as.numeric(canopy_variance),
  cumulative = cumsum(as.numeric(canopy_variance))
)
write_table(canopy_variance_table, "04_canopy_pca_variance.csv")

canopy_scores_output <- data.frame(plot = rownames(canopy_scores), canopy_scores,
                                   row.names = NULL)
write_table(canopy_scores_output, "05_canopy_pca_scores.csv")

canopy_loadings_output <- data.frame(species = rownames(canopy_loadings),
                                     canopy_loadings, row.names = NULL)
write_table(canopy_loadings_output, "06_canopy_pca_loadings.csv")

canopy_axis_names <- paste0("canopy_PC", 1:8)
canopy_score_matrix <- as.matrix(canopy_scores)
colnames(canopy_score_matrix) <- canopy_axis_names
missing_canopy <- setdiff(metadata$Plot, rownames(canopy_score_matrix))
if (length(missing_canopy) > 0L) {
  stop("No canopy record for plots: ", paste(missing_canopy, collapse = ", "))
}
metadata[, canopy_axis_names] <- canopy_score_matrix[metadata$Plot, , drop = FALSE]
plot_metadata <- metadata[!duplicated(metadata$Plot), , drop = FALSE]

# -----------------------------------------------------------------------------
# 6. Coverage-based diversity and the 0.90-0.99 coverage scan
# -----------------------------------------------------------------------------

coverage_levels <- seq(0.90, 0.99, by = 0.01)
target_coverage <- 0.99

phylum <- clean_phylum_names(tax[, "Phylum"])
sample_relative_abundance <- otu / rowSums(otu)
mean_relative_by_taxon <- colMeans(sample_relative_abundance)
mean_relative_by_phylum <- tapply(mean_relative_by_taxon, phylum, sum, na.rm = TRUE)
mean_relative_by_phylum <- sort(mean_relative_by_phylum, decreasing = TRUE)
major_phyla <- names(mean_relative_by_phylum)[seq_len(10)]

write_table(data.frame(
  rank = seq_along(mean_relative_by_phylum),
  phylum = names(mean_relative_by_phylum),
  mean_relative_abundance = as.numeric(mean_relative_by_phylum)
), "07_phylum_mean_relative_abundance.csv")

coverage_scan_rows <- list()
coverage_diagnostic_rows <- list()
for (habitat_name in c("litter", "soil")) {
  habitat_samples <- rownames(metadata)[metadata$habitat == habitat_name]
  for (community in c("All bacteria", major_phyla)) {
    taxon_keep <- if (community == "All bacteria") {
      rep(TRUE, ncol(otu))
    } else {
      !is.na(phylum) & phylum == community
    }
    estimates <- coverage_estimates(otu[habitat_samples, taxon_keep, drop = FALSE],
                                    levels = coverage_levels)
    for (level in coverage_levels) {
      at_level <- estimates[abs(estimates$target_coverage - level) < 1e-9 &
                              estimates$retained, , drop = FALSE]
      factor_at_level <- at_level$size_m / at_level$reads
      coverage_scan_rows[[paste(habitat_name, community, level)]] <- data.frame(
        habitat = habitat_name,
        community = community,
        target_coverage = level,
        libraries_in_habitat = length(habitat_samples),
        libraries_with_more_than_one_read = estimates$n_libraries_in_subset[1],
        libraries_retained = nrow(at_level),
        median_extrapolation_factor = median(factor_at_level),
        maximum_extrapolation_factor = max(factor_at_level),
        libraries_beyond_twice_observed_reads = sum(factor_at_level > 2)
      )
    }
    at_target <- estimates[abs(estimates$target_coverage - target_coverage) < 1e-9, ]
    coverage_diagnostic_rows[[paste(habitat_name, community)]] <- data.frame(
      habitat = habitat_name,
      community = community,
      otus = at_target$n_otus_in_subset[1],
      libraries_in_habitat = length(habitat_samples),
      libraries_retained_at_0.99 = sum(at_target$retained),
      median_reads = median(at_target$reads),
      minimum_reads = min(at_target$reads),
      median_observed_coverage = median(at_target$observed_coverage),
      minimum_observed_coverage = min(at_target$observed_coverage),
      median_extrapolation_factor_at_0.99 = median(at_target$size_m[at_target$retained] /
                                                    at_target$reads[at_target$retained]),
      libraries_beyond_twice_observed_reads_at_0.99 =
        sum(at_target$size_m[at_target$retained] / at_target$reads[at_target$retained] > 2)
    )
  }
}
coverage_scan <- do.call(rbind, coverage_scan_rows)
coverage_diagnostics <- do.call(rbind, coverage_diagnostic_rows)
write_table(coverage_scan, "08_coverage_scan_0.90_0.99.csv")
write_table(coverage_diagnostics, "09_coverage_diagnostics_at_0.99.csv")

# Highest coverage in the scan that retains every library in both habitats.
overall_scan <- coverage_scan[coverage_scan$community == "All bacteria", ]
levels_all_retained <- tapply(
  overall_scan$libraries_retained == overall_scan$libraries_in_habitat,
  overall_scan$target_coverage, all
)
stopifnot(isTRUE(levels_all_retained[[as.character(target_coverage)]]))

# Highest coverage that every library reaches within twice its observed read
# count, used for a sensitivity analysis of the four primary tests.
levels_within_double <- tapply(
  overall_scan$libraries_beyond_twice_observed_reads == 0 &
    overall_scan$libraries_retained == overall_scan$libraries_in_habitat,
  overall_scan$target_coverage, all
)
sensitivity_coverage <- if (any(levels_within_double)) {
  max(as.numeric(names(levels_within_double)[levels_within_double]))
} else {
  NA_real_
}

coverage_data <- list()
for (habitat_name in c("litter", "soil")) {
  sample_ids <- rownames(metadata)[metadata$habitat == habitat_name]
  coverage_data[[habitat_name]] <- estimate_coverage_q1(
    otu[sample_ids, , drop = FALSE], metadata, coverage = target_coverage
  )
  coverage_data[[habitat_name]]$habitat <- habitat_name
}

coverage_all <- do.call(rbind, coverage_data)
rownames(coverage_all) <- coverage_all$sample_id
write_table(coverage_all[, c("sample_id", "habitat", "q1", "sqrt_q1",
                            "reads_in_subset", "observed_coverage", "size_m",
                            "extrapolation_factor", "Elevation_relative",
                            "Elevation", "Latitude", "Location", "PermGroup",
                            "Mountain_range")],
            "10_coverage_q1_overall.csv")

# Habitat comparison repeated with coverage-standardized q = 1 diversity.
coverage_habitat_model <- lmer(sqrt_q1 ~ habitat + (1 | Location),
                               data = coverage_all, REML = TRUE)
coverage_habitat_anova <- anova(coverage_habitat_model)
coverage_habitat_table <- data.frame(
  response = "sqrt(q = 1) at 0.99 coverage",
  n_litter = sum(coverage_all$habitat == "litter"),
  n_soil = sum(coverage_all$habitat == "soil"),
  mean_q1_litter = mean(coverage_all$q1[coverage_all$habitat == "litter"]),
  mean_q1_soil = mean(coverage_all$q1[coverage_all$habitat == "soil"]),
  soil_minus_litter_estimate = fixef(coverage_habitat_model)[["habitatsoil"]],
  F_value = coverage_habitat_anova[["F value"]][1],
  denominator_df = coverage_habitat_anova[["DenDF"]][1],
  p = coverage_habitat_anova[["Pr(>F)"]][1]
)
write_table(coverage_habitat_table, "11_habitat_comparison_coverage_0.99.csv")

# -----------------------------------------------------------------------------
# 7. Figure 2: geographic models, BH correction, raw-elevation sensitivity
# -----------------------------------------------------------------------------

elevation_variables <- c(relative = "Elevation_relative", raw = "Elevation")
geographic_fits <- list()
geographic_model_rows <- list()
geographic_aic_rows <- list()
for (habitat_name in names(coverage_data)) {
  dat <- coverage_data[[habitat_name]]
  for (elevation_label in names(elevation_variables)) {
    elevation_variable <- unname(elevation_variables[[elevation_label]])
    for (with_interaction in c(FALSE, TRUE)) {
      key <- paste(habitat_name, elevation_label, with_interaction, sep = "_")
      fit <- fit_candidate_lmms(
        dat, response = "sqrt_q1", elevation_variable = elevation_variable,
        interaction = with_interaction
      )
      geographic_fits[[key]] <- fit
      geographic_model_rows[[key]] <- model_summary_row(
        fit, habitat_name, elevation_label, with_interaction, nrow(dat)
      )
      geographic_aic_rows[[key]] <- candidate_aic_table(
        fit, habitat_name, elevation_label, with_interaction
      )
    }
  }
}

geographic_model_table <- do.call(rbind, geographic_model_rows)
geographic_aic_table <- do.call(rbind, geographic_aic_rows)
# Overall AIC ranking of the eight candidates (additive and interaction).
geographic_aic_table$delta_AIC_all_eight <- ave(
  geographic_aic_table$AIC,
  paste(geographic_aic_table$habitat, geographic_aic_table$elevation_type),
  FUN = function(x) x - min(x)
)
geographic_aic_table$selected_among_all_eight <- geographic_aic_table$delta_AIC_all_eight == 0

# Benjamini-Hochberg correction across the four primary habitat x gradient
# tests (relative elevation, additive models), and separately across the four
# raw-elevation sensitivity tests.
geographic_model_table$elevation_p_BH <- NA_real_
geographic_model_table$latitude_p_BH <- NA_real_
for (elevation_label in names(elevation_variables)) {
  rows <- geographic_model_table$elevation_type == elevation_label &
    !geographic_model_table$interaction
  adjusted <- p.adjust(c(geographic_model_table$elevation_p[rows],
                         geographic_model_table$latitude_p[rows]), method = "BH")
  geographic_model_table$elevation_p_BH[rows] <- adjusted[seq_len(sum(rows))]
  geographic_model_table$latitude_p_BH[rows] <- adjusted[sum(rows) + seq_len(sum(rows))]
}
write_table(geographic_model_table, "12_geographic_LMM_results.csv")
write_table(geographic_aic_table, "13_geographic_LMM_AIC_selection.csv")

primary_rows <- geographic_model_table$elevation_type == "relative" &
  !geographic_model_table$interaction
stopifnot(all(geographic_aic_table$model[geographic_aic_table$selected_among_all_eight &
                                           geographic_aic_table$elevation_type == "relative"] ==
                geographic_model_table$selected_model[primary_rows]))

relative_additive_fits <- list(
  litter = geographic_fits[["litter_relative_FALSE"]],
  soil = geographic_fits[["soil_relative_FALSE"]]
)
primary_p_BH <- list(
  litter = geographic_model_table[primary_rows & geographic_model_table$habitat == "litter", ],
  soil = geographic_model_table[primary_rows & geographic_model_table$habitat == "soil", ]
)

# Sensitivity: the four primary tests at the highest coverage reached by every
# library within twice its observed read count.
sensitivity_rows <- list()
if (!is.na(sensitivity_coverage)) {
  for (habitat_name in c("litter", "soil")) {
    sample_ids <- rownames(metadata)[metadata$habitat == habitat_name]
    dat <- estimate_coverage_q1(otu[sample_ids, , drop = FALSE], metadata,
                                coverage = sensitivity_coverage)
    fit <- fit_candidate_lmms(dat, response = "sqrt_q1",
                              elevation_variable = "Elevation_relative")
    sensitivity_rows[[habitat_name]] <- model_summary_row(
      fit, habitat_name, "relative", FALSE, nrow(dat)
    )
    sensitivity_rows[[habitat_name]]$target_coverage <- sensitivity_coverage
    sensitivity_rows[[habitat_name]]$correlation_with_0.99_estimates <- cor(
      dat$q1, coverage_data[[habitat_name]][dat$sample_id, "q1"]
    )
  }
  coverage_sensitivity_table <- do.call(rbind, sensitivity_rows)
  adjusted <- p.adjust(c(coverage_sensitivity_table$elevation_p,
                         coverage_sensitivity_table$latitude_p), method = "BH")
  coverage_sensitivity_table$elevation_p_BH <- adjusted[1:2]
  coverage_sensitivity_table$latitude_p_BH <- adjusted[3:4]
  write_table(coverage_sensitivity_table, "14_geographic_LMM_coverage_sensitivity.csv")
}

# Colours of the mountain ranges (colour-blind-friendly palette).
mountain_range_colors <- c(
  Alps = "#000000", Baetic_System = "#E69F00", Carpathian_Mountains = "#56B4E9",
  Dinaric_Alps = "#009E73", Harz_Mountains = "#F0E442",
  Iceland_Mountains = "#0072B2", Pyrenees = "#D55E00", Rila_Rhodopes = "#CC79A7",
  Scandinavian_Mountains = "#999999", Taygetos_Mountains = "#44AA99"
)
stopifnot(all(unique(metadata$Mountain_range) %in% names(mountain_range_colors)))

make_geographic_panel <- function(data, gradient, fit, p_adjusted) {
  label <- NULL
  significant <- !is.na(p_adjusted) && p_adjusted < 0.05
  if (significant) {
    label <- sprintf(
      "R²m = %.3f\nR²c = %.3f\n%s",
      fit$r2[["marginal"]], fit$r2[["conditional"]],
      format_p(p_adjusted, "p(BH)")
    )
  }

  panel <- ggplot(data, aes(x = .data[[gradient]], y = q1,
                            color = Mountain_range))
  if (significant) {
    trend <- lmm_trend(fit$best_model, data, gradient)
    panel <- panel +
      geom_ribbon(data = trend, aes(x = x, ymin = lower, ymax = upper),
                  inherit.aes = FALSE, fill = "grey60", alpha = 0.35) +
      geom_line(data = trend, aes(x = x, y = q1), inherit.aes = FALSE,
                color = "black", linewidth = 0.5)
  }
  panel <- panel +
    geom_point(alpha = 0.7, size = 2.4, shape = 16) +
    scale_color_manual(values = mountain_range_colors,
                       labels = function(x) gsub("_", " ", x),
                       limits = names(mountain_range_colors)) +
    labs(
      x = if (gradient == "Elevation_relative") "Elevation" else "Latitude",
      y = "Hill number q=1",
      color = "Mountain Range"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8))) +
    theme_bw(base_size = 11, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.tag = element_text(size = 18),
      plot.tag.position = c(0.01, 0.985),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 15)
    )
  if (significant) {
    panel <- panel +
      annotate(
        "text",
        x = if (gradient == "Elevation_relative") min(data[[gradient]]) else max(data[[gradient]]),
        y = max(data$q1), label = label, family = "serif",
        hjust = if (gradient == "Elevation_relative") 0 else 1, vjust = 1, size = 4.5
      )
  }
  panel
}

fig2a <- make_geographic_panel(coverage_data$litter, "Elevation_relative",
                               relative_additive_fits$litter,
                               primary_p_BH$litter$elevation_p_BH)
fig2b <- make_geographic_panel(coverage_data$litter, "Latitude",
                               relative_additive_fits$litter,
                               primary_p_BH$litter$latitude_p_BH)
fig2c <- make_geographic_panel(coverage_data$soil, "Elevation_relative",
                               relative_additive_fits$soil,
                               primary_p_BH$soil$elevation_p_BH)
fig2d <- make_geographic_panel(coverage_data$soil, "Latitude",
                               relative_additive_fits$soil,
                               primary_p_BH$soil$latitude_p_BH)

# Rows: litter (a, b) and soil (c, d); columns: elevation and latitude.
habitat_label <- function(text) {
  wrap_elements(full = grid::textGrob(
    text, rot = 90, gp = grid::gpar(fontfamily = "serif", fontsize = 20)
  ))
}
fig2 <- wrap_plots(
  habitat_label("Litter"),
  fig2a + labs(tag = "a)", title = "Elevation"),
  fig2b + labs(tag = "b)", title = "Latitude"),
  habitat_label("Soil"),
  fig2c + labs(tag = "c)", title = " "),
  fig2d + labs(tag = "d)", title = " "),
  ncol = 3, widths = c(0.06, 1, 1), guides = "collect"
) & theme(legend.position = "right")
ggsave(file.path(figure_dir, "Fig2.pdf"), fig2,
       width = 14.5, height = 8.8, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig2.png"), fig2,
       width = 14.5, height = 8.8, dpi = 180)

# -----------------------------------------------------------------------------
# 8. Residual analyses: climate, soil, and combined environmental models
# -----------------------------------------------------------------------------

climate_variables <- c(
  "Temperature_Seasonality", "Min_Temp_Coldest_Month",
  "Mean_Temp_Wettest_Quarter", "Mean_Temp_Driest_Quarter",
  "Precipitation_Seasonality", "Precipitation_Warmest_Quarter"
)
soil_variables <- c("n", "Ca", "pH", "p_MEL")
environment_sets <- list(
  climate = climate_variables,
  soil = soil_variables,
  combined = c(climate_variables, soil_variables)
)

run_residual_analysis <- function(data, response = "sqrt_q1",
                                  elevation_variable = "Elevation_relative") {
  rows <- list()
  models <- list()
  for (set_name in names(environment_sets)) {
    variables <- environment_sets[[set_name]]
    complete <- complete.cases(data[, c(response, variables, elevation_variable,
                                         "Latitude", "PermGroup", "Location")])
    model_data <- data[complete, , drop = FALSE]
    environment_formula <- as.formula(sprintf(
      "%s ~ %s", response, paste(variables, collapse = " + ")
    ))
    environment_model <- lm(environment_formula, data = model_data)
    model_data$environment_residual <- residuals(environment_model)
    geographic_fit <- fit_candidate_lmms(
      model_data, response = "environment_residual",
      elevation_variable = elevation_variable, interaction = FALSE
    )
    rows[[set_name]] <- cbind(
      data.frame(
        environmental_set = set_name,
        n_samples = nrow(model_data),
        environmental_R2 = summary(environment_model)$r.squared,
        selected_geographic_model = geographic_fit$best_name,
        singular_fit = geographic_fit$singular
      ),
      gradient_columns(geographic_fit$elevation, "residual_elevation"),
      gradient_columns(geographic_fit$latitude, "residual_latitude")
    )
    models[[set_name]] <- list(
      environment_model = environment_model,
      geographic_fit = geographic_fit
    )
  }
  list(table = do.call(rbind, rows), models = models)
}

overall_residual_results <- list()
overall_residual_rows <- list()
for (habitat_name in names(coverage_data)) {
  for (elevation_label in names(elevation_variables)) {
    result <- run_residual_analysis(
      coverage_data[[habitat_name]],
      elevation_variable = unname(elevation_variables[[elevation_label]])
    )
    key <- paste(habitat_name, elevation_label, sep = "_")
    result$table <- cbind(habitat = habitat_name, elevation_type = elevation_label,
                          result$table)
    overall_residual_results[[key]] <- result
    overall_residual_rows[[key]] <- result$table
  }
}
overall_residual_table <- do.call(rbind, overall_residual_rows)
write_table(overall_residual_table, "15_overall_residual_analysis.csv")

# Environmental correlation matrix for Fig. S1. The figure shows the complete
# candidate set before collinearity filtering, ordered as climate, canopy, and
# soil. Environmental measurements occur twice in sample_data (paired litter
# and soil); use one row per plot here.
rf_climate_candidates <- c(
  "Mean_Diurnal_Range", "MAT", "MAP",
  "Temperature_Seasonality", "Max_Temp_Warmest_Month", "Min_Temp_Coldest_Month",
  "Temperature_Annual_Range", "Mean_Temp_Wettest_Quarter",
  "Mean_Temp_Driest_Quarter", "Mean_Temp_Warmest_Quarter",
  "Mean_Temp_Coldest_Quarter", "Precipitation_Wettest_Month",
  "Precipitation_Driest_Month", "Precipitation_Seasonality",
  "Precipitation_Wettest_Quarter", "Precipitation_Driest_Quarter",
  "Precipitation_Warmest_Quarter", "Precipitation_Coldest_Quarter"
)
rf_soil_candidates <- c("c", "n", "Ca", "pH", "p_MEL")
rf_canopy_candidates <- canopy_axis_names
correlation_variables <- c(
  rf_climate_candidates, rf_canopy_candidates, rf_soil_candidates
)
correlation_data <- plot_metadata[, correlation_variables, drop = FALSE]
correlation_matrix <- cor(correlation_data, use = "pairwise.complete.obs")
write_table(data.frame(variable = rownames(correlation_matrix), correlation_matrix,
                       check.names = FALSE),
            "16_environmental_correlation_matrix.csv")
ca_ph_correlation <- correlation_matrix["Ca", "pH"]

correlation_labels <- c(
  setNames(paste0("PC", 1:8, " canopy"), canopy_axis_names),
  c = "C", n = "N", p_MEL = "P"
)
correlation_plot_matrix <- correlation_matrix
pretty_names <- ifelse(rownames(correlation_plot_matrix) %in% names(correlation_labels),
                       correlation_labels[rownames(correlation_plot_matrix)],
                       gsub("_", " ", rownames(correlation_plot_matrix)))
dimnames(correlation_plot_matrix) <- list(pretty_names, pretty_names)

grDevices::pdf(file.path(figure_dir, "FigS1_correlation.pdf"),
               width = 11, height = 9)
corrplot::corrplot(correlation_plot_matrix, method = "circle", type = "upper",
                   tl.col = "black", tl.cex = 0.7)
grDevices::dev.off()

grDevices::png(file.path(figure_dir, "FigS1_correlation.png"),
               width = 11, height = 9, units = "in", res = 220)
corrplot::corrplot(correlation_plot_matrix, method = "circle", type = "upper",
                   tl.col = "black", tl.cex = 0.7)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# 9. Phylum-level coverage, model selection, residual analyses, Figures 3-4
# -----------------------------------------------------------------------------
# Eligibility is tested for every phylum-habitat combination with the same
# rule: at least two OTUs and finite q = 1 estimates at 0.99 coverage for at
# least 90% of the libraries of that habitat.

minimum_retained_fraction <- 0.90
phylum_analysis <- list()
phylum_rows <- list()

for (habitat_name in c("soil", "litter")) {
  habitat_samples <- rownames(metadata)[metadata$habitat == habitat_name]
  for (phylum_name in major_phyla) {
    taxon_keep <- !is.na(phylum) & phylum == phylum_name
    phylum_data <- estimate_coverage_q1(
      otu[habitat_samples, taxon_keep, drop = FALSE], metadata,
      coverage = target_coverage
    )
    eligible <- sum(taxon_keep) >= 2L &&
      nrow(phylum_data) >= minimum_retained_fraction * length(habitat_samples)
    key <- paste(phylum_name, habitat_name, sep = "_")
    if (!eligible) {
      phylum_rows[[key]] <- data.frame(phylum = phylum_name, habitat = habitat_name,
                                       n_samples = nrow(phylum_data), eligible = FALSE)
      next
    }
    phylum_data$habitat <- habitat_name
    fit <- fit_candidate_lmms(phylum_data, response = "sqrt_q1",
                              elevation_variable = "Elevation_relative",
                              interaction = FALSE)
    phylum_analysis[[key]] <- list(data = phylum_data, fit = fit)
    phylum_rows[[key]] <- cbind(
      data.frame(
        phylum = phylum_name,
        habitat = habitat_name,
        n_samples = nrow(phylum_data),
        eligible = TRUE,
        median_extrapolation_factor = median(phylum_data$extrapolation_factor),
        coverage = target_coverage,
        selected_model = fit$best_name,
        AIC = min(fit$aic),
        singular_fit = fit$singular
      ),
      gradient_columns(fit$elevation, "elevation"),
      gradient_columns(fit$latitude, "latitude"),
      data.frame(
        R2_marginal = fit$r2[["marginal"]],
        R2_conditional = fit$r2[["conditional"]]
      )
    )
  }
}

phylum_results <- do.call(rbind, lapply(phylum_rows, function(x) {
  template <- phylum_rows[[which(vapply(phylum_rows, ncol, integer(1)) ==
                                   max(vapply(phylum_rows, ncol, integer(1))))[1]]]
  missing <- setdiff(names(template), names(x))
  x[missing] <- NA
  x[, names(template)]
}))
phylum_results$elevation_p_BH <- p.adjust(phylum_results$elevation_p, method = "BH")
phylum_results$latitude_p_BH <- p.adjust(phylum_results$latitude_p, method = "BH")
# Curvature: quadratic coefficients of the selected polynomial models, BH
# adjusted across all quadratic terms retained by model selection.
quadratic_p <- c(phylum_results$elevation_quadratic_p, phylum_results$latitude_quadratic_p)
quadratic_p_BH <- rep(NA_real_, length(quadratic_p))
quadratic_p_BH[!is.na(quadratic_p)] <- p.adjust(quadratic_p[!is.na(quadratic_p)],
                                                method = "BH")
phylum_results$elevation_quadratic_p_BH <- quadratic_p_BH[seq_len(nrow(phylum_results))]
phylum_results$latitude_quadratic_p_BH <- quadratic_p_BH[nrow(phylum_results) +
                                                           seq_len(nrow(phylum_results))]
write_table(phylum_results, "17_phylum_LMM_results.csv")

# Residual analyses are required only for phylum-habitat combinations with at
# least one significant geographic relationship in the main LMM analysis.
significant_phylum_rows <- phylum_results[
  (!is.na(phylum_results$elevation_p_BH) & phylum_results$elevation_p_BH < 0.05) |
  (!is.na(phylum_results$latitude_p_BH) & phylum_results$latitude_p_BH < 0.05),
]

phylum_residual_rows <- list()
for (i in seq_len(nrow(significant_phylum_rows))) {
  row <- significant_phylum_rows[i, ]
  key <- paste(row$phylum, row$habitat, sep = "_")
  residual_result <- run_residual_analysis(phylum_analysis[[key]]$data)
  residual_result$table$phylum <- row$phylum
  residual_result$table$habitat <- row$habitat
  residual_result$table$significant_gradients <- paste(c(
    if (!is.na(row$elevation_p_BH) && row$elevation_p_BH < 0.05) "elevation",
    if (!is.na(row$latitude_p_BH) && row$latitude_p_BH < 0.05) "latitude"
  ), collapse = "+")
  phylum_residual_rows[[key]] <- residual_result$table
}
phylum_residual_table <- do.call(rbind, phylum_residual_rows)
phylum_residual_table <- phylum_residual_table[
  , c("phylum", "habitat", "significant_gradients",
      setdiff(names(phylum_residual_table),
              c("phylum", "habitat", "significant_gradients")))
]
write_table(phylum_residual_table, "18_phylum_residual_analysis.csv")

make_phylum_figure <- function(gradient = c("Elevation_relative", "Latitude")) {
  gradient <- match.arg(gradient)
  prefix <- if (gradient == "Elevation_relative") "elevation" else "latitude"
  p_column <- paste0(prefix, "_p_BH")
  significant_phyla <- unique(phylum_results$phylum[
    !is.na(phylum_results[[p_column]]) & phylum_results[[p_column]] < 0.05
  ])
  plot_rows <- list()
  annotation_rows <- list()
  smooth_layers <- list()

  for (key in names(phylum_analysis)) {
    entry <- phylum_analysis[[key]]
    phylum_name <- sub("_(soil|litter)$", "", key)
    habitat_name <- sub("^.*_", "", key)
    if (!phylum_name %in% significant_phyla) next
    result_row <- phylum_results[
      phylum_results$phylum == phylum_name & phylum_results$habitat == habitat_name,
    ]
    dat <- entry$data
    dat$phylum <- phylum_name
    dat$habitat <- habitat_name
    significant <- !is.na(result_row[[p_column]]) && result_row[[p_column]] < 0.05
    dat$significant <- significant
    plot_rows[[key]] <- dat[, c(gradient, "q1", "phylum", "habitat", "significant")]

    if (significant) {
      annotation_rows[[key]] <- data.frame(
        phylum = phylum_name,
        habitat = habitat_name,
        x = if (habitat_name == "soil") {
          min(dat[[gradient]]) + 0.02 * diff(range(dat[[gradient]]))
        } else {
          max(dat[[gradient]]) - 0.02 * diff(range(dat[[gradient]]))
        },
        y = max(dat$q1) - 0.02 * diff(range(dat$q1)),
        hjust = if (habitat_name == "soil") 0 else 1,
        label = sprintf(
          "R²m = %.3f\nR²c = %.3f\n%s",
          result_row$R2_marginal, result_row$R2_conditional,
          format_p(result_row[[p_column]], "p(BH)")
        )
      )
      trend <- lmm_trend(entry$fit$best_model, dat, gradient)
      trend$phylum <- phylum_name
      trend$habitat <- factor(habitat_name, levels = c("litter", "soil"))
      smooth_layers[[paste0(key, "_band")]] <- geom_ribbon(
        data = trend, aes(x = x, ymin = lower, ymax = upper, fill = habitat),
        inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE
      )
      smooth_layers[[paste0(key, "_line")]] <- geom_line(
        data = trend, aes(x = x, y = q1, color = habitat),
        inherit.aes = FALSE, linewidth = 0.8, show.legend = FALSE
      )
    }
  }

  plot_data <- do.call(rbind, plot_rows)
  annotations <- do.call(rbind, annotation_rows)
  plot_data$habitat <- factor(plot_data$habitat, levels = c("litter", "soil"))
  # Panels in alphabetical order, independent of the system locale.
  phylum_levels <- sort(unique(plot_data$phylum), method = "radix")
  plot_data$phylum <- factor(plot_data$phylum, levels = phylum_levels)
  if (!is.null(annotations)) {
    annotations$phylum <- factor(annotations$phylum, levels = phylum_levels)
  }

  figure <- ggplot(plot_data, aes(x = .data[[gradient]], y = q1,
                                  color = habitat)) +
    geom_point(aes(alpha = significant), size = 1.5) +
    scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.25), guide = "none") +
    scale_color_manual(values = c(litter = "#01665e", soil = "#8c510a"),
                       labels = c(litter = "Litter", soil = "Soil")) +
    scale_fill_manual(values = c(litter = "#01665e", soil = "#8c510a"), guide = "none") +
    facet_wrap(~phylum, scales = "free_y", ncol = 3) +
    labs(
      x = if (gradient == "Elevation_relative") "Elevation" else "Latitude (°)",
      y = "Hill number q=1", color = "Sample Type"
    ) +
    theme_figures(10) +
    theme(legend.position = "bottom")
  for (layer in smooth_layers) {
    layer$data$phylum <- factor(layer$data$phylum, levels = phylum_levels)
    figure <- figure + layer
  }

  if (!is.null(annotations) && nrow(annotations) > 0L) {
    figure <- figure + geom_label(
      data = annotations,
      aes(x = x, y = y, label = label, color = habitat, hjust = hjust),
      vjust = 1, fill = "white", linewidth = 0, size = 3,
      inherit.aes = FALSE, show.legend = FALSE
    )
  }
  figure
}

fig3 <- make_phylum_figure("Elevation_relative")
fig4 <- make_phylum_figure("Latitude")
ggsave(file.path(figure_dir, "Fig3.pdf"), fig3,
       width = 13, height = 9, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig3.png"), fig3,
       width = 13, height = 9, dpi = 180)
ggsave(file.path(figure_dir, "Fig4.pdf"), fig4,
       width = 13, height = 9, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig4.png"), fig4,
       width = 13, height = 9, dpi = 180)

# -----------------------------------------------------------------------------
# 10. Random Forest: tuning, OOB performance, importance, PDPs, Figure 5
# -----------------------------------------------------------------------------

find_correlated_variables <- function(correlation_matrix, threshold = 0.7) {
  if (ncol(correlation_matrix) < 2L) return(character())
  remove <- character()
  for (i in seq_len(ncol(correlation_matrix) - 1L)) {
    for (j in (i + 1L):ncol(correlation_matrix)) {
      if (abs(correlation_matrix[i, j]) > threshold) {
        mean_i <- mean(abs(correlation_matrix[i, -i]))
        mean_j <- mean(abs(correlation_matrix[j, -j]))
        remove <- c(remove, if (mean_i > mean_j) {
          colnames(correlation_matrix)[i]
        } else {
          colnames(correlation_matrix)[j]
        })
      }
    }
  }
  unique(remove)
}

reduce_collinearity <- function(data, variables, threshold = 0.7) {
  variables <- intersect(variables, names(data))
  correlation_matrix <- cor(data[, variables, drop = FALSE],
                            use = "pairwise.complete.obs")
  setdiff(variables, find_correlated_variables(correlation_matrix, threshold))
}

rf_data <- metadata
# Response: Shannon diversity of the rarefied table used for Fig. 1.
rf_data$shannon_rarefied <- shannon_rarefied[rownames(rf_data)]

rf_climate_selected <- reduce_collinearity(rf_data, rf_climate_candidates)
rf_soil_selected <- reduce_collinearity(rf_data, rf_soil_candidates)
rf_canopy_selected <- reduce_collinearity(rf_data, rf_canopy_candidates)
stopifnot(setequal(rf_climate_selected, climate_variables),
          setequal(rf_soil_selected, soil_variables))
rf_predictors <- c(rf_climate_selected, rf_soil_selected, rf_canopy_selected)

predictor_category <- function(x) {
  ifelse(x %in% rf_climate_selected, "Climate",
         ifelse(x %in% rf_soil_selected, "Soil", "Canopy"))
}

rf_predictor_table <- data.frame(
  variable = rf_predictors,
  category = predictor_category(rf_predictors)
)
write_table(rf_predictor_table, "19_random_forest_selected_predictors.csv")

tune_and_fit_rf <- function(model_data, response, predictors, ntree_final = 1000,
                            importance = TRUE) {
  formula <- as.formula(paste(response, "~ ."))
  model_data <- model_data[, c(response, predictors), drop = FALSE]
  mtry_values <- seq(
    floor(sqrt(length(predictors))),
    min(ceiling(length(predictors) / 3), length(predictors)), by = 1
  )
  tuning <- do.call(rbind, lapply(mtry_values, function(mtry_value) {
    set.seed(123)
    tuning_model <- randomForest(formula, data = model_data, mtry = mtry_value,
                                 ntree = 500, importance = FALSE)
    data.frame(mtry = mtry_value, OOB_MSE = tail(tuning_model$mse, 1))
  }))
  best_mtry <- tuning$mtry[which.min(tuning$OOB_MSE)]
  set.seed(456)
  model <- randomForest(formula, data = model_data, mtry = best_mtry,
                        ntree = ntree_final, importance = importance)
  list(model = model, model_data = model_data, tuning = tuning)
}

run_random_forest <- function(data, habitat_name, response, predictors) {
  model_data <- data[complete.cases(data[, c(response, predictors)]), , drop = FALSE]
  fitted <- tune_and_fit_rf(model_data, response, predictors)
  model <- fitted$model
  importance_matrix <- randomForest::importance(model, type = 1)
  importance <- data.frame(
    habitat = habitat_name,
    variable = rownames(importance_matrix),
    importance = importance_matrix[, 1],
    category = predictor_category(rownames(importance_matrix))
  )
  importance <- importance[order(importance$importance, decreasing = TRUE), ]
  importance$rank <- seq_len(nrow(importance))

  oob_mse <- tail(model$mse, 1)
  top10 <- head(importance, 10)
  performance <- data.frame(
    habitat = habitat_name,
    response = response,
    n_samples = nrow(model_data),
    n_predictors = length(predictors),
    ntree = model$ntree,
    mtry = model$mtry,
    OOB_RMSE = sqrt(oob_mse),
    OOB_R2 = 1 - oob_mse / mean((model$y - mean(model$y))^2),
    top10_canopy = sum(top10$category == "Canopy"),
    top10_climate = sum(top10$category == "Climate"),
    top10_soil = sum(top10$category == "Soil")
  )

  list(model = model, model_data = fitted$model_data,
       tuning = cbind(habitat = habitat_name, fitted$tuning),
       importance = importance, performance = performance)
}

rf_results <- list()
for (habitat_name in c("litter", "soil")) {
  rf_results[[habitat_name]] <- run_random_forest(
    rf_data[rf_data$habitat == habitat_name, , drop = FALSE], habitat_name,
    "shannon_rarefied", rf_predictors
  )
}

flatten_rf <- function(field) {
  do.call(rbind, lapply(rf_results, `[[`, field))
}
rf_tuning <- flatten_rf("tuning")
rf_importance <- flatten_rf("importance")
rf_performance <- flatten_rf("performance")
write_table(rf_tuning, "20_random_forest_tuning.csv")
write_table(rf_performance, "21_random_forest_performance.csv")
write_table(rf_importance, "22_random_forest_importance.csv")

variable_labels <- c(
  setNames(paste0("PC", 1:8, " canopy"), canopy_axis_names),
  p_MEL = "P", Min_Temp_Coldest_Month = "Min Temp Coldest Month",
  Temperature_Seasonality = "Temperature Seasonality",
  Mean_Temp_Wettest_Quarter = "Mean Temp Wettest Quarter",
  Mean_Temp_Driest_Quarter = "Mean Temp Driest Quarter",
  Precipitation_Warmest_Quarter = "Precipitation Warmest Quarter",
  Precipitation_Seasonality = "Precipitation Seasonality", n = "N"
)
label_variable <- function(x) {
  ifelse(x %in% names(variable_labels), variable_labels[x], x)
}

rf_top10 <- do.call(rbind, lapply(split(rf_importance, rf_importance$habitat),
                                  function(x) head(x, 10)))
rf_top10$label <- label_variable(rf_top10$variable)
rf_top10$plot_label <- paste(rf_top10$habitat, rf_top10$label, sep = "___")
rf_top10$plot_label <- factor(
  rf_top10$plot_label,
  levels = rf_top10$plot_label[order(rf_top10$habitat, rf_top10$importance)]
)

fig5 <- ggplot(rf_top10, aes(plot_label, importance, fill = category)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~habitat, scales = "free_y",
             labeller = as_labeller(c(litter = "Litter", soil = "Soil"))) +
  scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(values = c(Canopy = "#228B22", Climate = "#4682B4",
                               Soil = "#8B4513")) +
  labs(x = "Environmental Variables", y = "Variable Importance (%IncMSE)",
       fill = "Category") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

ggsave(file.path(figure_dir, "Fig5.pdf"), fig5,
       width = 11, height = 8, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig5.png"), fig5,
       width = 11, height = 8, dpi = 200)

if (run_pdp) {
  make_pdp_plot <- function(result, variable) {
    partial_result <- pdp::partial(
      result$model, pred.var = variable, train = result$model_data,
      progress = "none"
    )
    autoplot(partial_result, rug = TRUE, train = result$model_data) +
      labs(x = label_variable(variable), y = "Partial effect") +
      theme_minimal(base_size = 9)
  }

  for (habitat_name in c("litter", "soil")) {
    result <- rf_results[[habitat_name]]
    variables <- head(result$importance$variable, 10)
    plots <- lapply(variables, function(variable) make_pdp_plot(result, variable))
    pdp_figure <- wrap_plots(plots, ncol = 2) +
      plot_annotation(title = paste("Partial dependence -", tools::toTitleCase(habitat_name)))
    figure_number <- if (habitat_name == "litter") "S2" else "S3"
    ggsave(file.path(figure_dir,
                     paste0("Fig", figure_number, "_PDP_", habitat_name, ".pdf")),
           pdp_figure, width = 8, height = 10, device = grDevices::pdf)
    ggsave(file.path(figure_dir,
                     paste0("Fig", figure_number, "_PDP_", habitat_name, ".png")),
           pdp_figure, width = 8, height = 10, dpi = 220)
  }
}

# Species loadings for canopy axes appearing among either habitat's top ten.
important_canopy_axes <- unique(sub("^canopy_", "", rf_top10$variable[
  rf_top10$category == "Canopy"
]))
important_canopy_axes <- important_canopy_axes[order(as.integer(sub("PC", "", important_canopy_axes)))]
loading_rows <- list()
for (axis_name in important_canopy_axes) {
  axis_values <- canopy_loadings[[axis_name]]
  top <- order(abs(axis_values), decreasing = TRUE)[seq_len(min(10, length(axis_values)))]
  loading_rows[[axis_name]] <- data.frame(
    axis = axis_name,
    species = rownames(canopy_loadings)[top],
    loading = axis_values[top],
    absolute_loading = abs(axis_values[top])
  )
}
important_loadings <- do.call(rbind, loading_rows)

# Species classes and colors follow Fig. S4. They describe tree
# type, not the sign of the PCA loading.
coniferous_species <- c(
  "Picea abies", "Pinus sylvestris", "Pinus sp.", "Abies alba",
  "Pinaceae sp.", "Pinus nigra", "Larix decidua", "Pinus cembra"
)
deciduous_species <- c(
  "Betula pendula", "Betula pubescens", "Sorbus aucuparia",
  "Corylus avellana", "Fraxinus excelsior", "Alnus glutinosa",
  "Alnus incana", "Salix sp.", "Fagus sylvatica", "Rhamnus cathartica",
  "Acer sp.", "Quercus sp.", "Aesculus hippocastanum", "Carpinus betulus",
  "Quercus rotundifolia", "Sorbus domestica", "Quercus faginea",
  "Populus tremula"
)
important_loadings$species_type <- ifelse(
  important_loadings$species %in% coniferous_species, "Coniferous",
  ifelse(important_loadings$species %in% deciduous_species,
         "Deciduous", "Other")
)
write_table(important_loadings, "23_canopy_top_species_loadings.csv")

loading_plots <- lapply(important_canopy_axes, function(axis_name) {
  x <- important_loadings[important_loadings$axis == axis_name, ]
  x$species <- factor(x$species, levels = x$species[order(x$absolute_loading)])
  axis_index <- match(axis_name, names(canopy_variance))
  ggplot(x, aes(species, loading, fill = species_type)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(
      values = c(Coniferous = "#228B22", Deciduous = "#D2691E",
                 Other = "grey50"), guide = "none"
    ) +
    labs(
      title = sprintf("Canopy %s (%.1f%%)", axis_name,
                      100 * canopy_variance[axis_index]),
      x = NULL, y = "Species contribution (PCA loading)"
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(face = "italic"))
})
fig_s4 <- wrap_plots(loading_plots, ncol = 2)
ggsave(file.path(figure_dir, "FigS4_canopy_loadings.pdf"), fig_s4,
       width = 10, height = max(6, ceiling(length(loading_plots) / 2) * 3.5),
       device = grDevices::pdf)
ggsave(file.path(figure_dir, "FigS4_canopy_loadings.png"), fig_s4,
       width = 10, height = max(6, ceiling(length(loading_plots) / 2) * 3.5),
       dpi = 220)

# -----------------------------------------------------------------------------
# 11. Save fitted models, software versions, and run summary
# -----------------------------------------------------------------------------

saveRDS(
  list(
    habitat_model = habitat_model,
    coverage_habitat_model = coverage_habitat_model,
    geographic_fits = geographic_fits,
    residual_results = overall_residual_results,
    phylum_analysis = phylum_analysis,
    rf_results = rf_results,
    canopy_pca = canopy_pca
  ),
  file.path(model_dir, "fitted_models.rds")
)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
write.csv(data.frame(
  package = c("R", unique(required_packages)),
  version = c(paste(R.version$major, R.version$minor, sep = "."),
              vapply(unique(required_packages), function(x) {
                as.character(utils::packageVersion(x))
              }, character(1)))
), file.path(output_dir, "package_versions.csv"), row.names = FALSE)

input_files <- c(phyloseq_file, canopy_file, elevation_reference_file)
write.csv(data.frame(
  file = basename(input_files),
  md5 = unname(tools::md5sum(input_files))
), file.path(output_dir, "input_checksums.csv"), row.names = FALSE)

geo_value <- function(habitat_name, elevation_label, field, interaction = FALSE) {
  geographic_model_table[
    geographic_model_table$habitat == habitat_name &
      geographic_model_table$elevation_type == elevation_label &
      geographic_model_table$interaction == interaction,
    field
  ][1]
}
residual_value <- function(habitat_name, elevation_label, set_name, field) {
  overall_residual_table[
    overall_residual_table$habitat == habitat_name &
      overall_residual_table$elevation_type == elevation_label &
      overall_residual_table$environmental_set == set_name, field
  ][1]
}

report_lines <- c(
  "RUN SUMMARY",
  "===========",
  paste("Run date:", as.character(Sys.time())),
  paste("R:", R.version.string),
  "",
  "Dataset:",
  sprintf("- %d samples (%d litter, %d soil)", nrow(otu),
          sum(metadata$habitat == "litter"), sum(metadata$habitat == "soil")),
  sprintf("- %s retained reads; depth %s-%s", format(sum(otu), big.mark = ","),
          format(min(rowSums(otu)), big.mark = ","),
          format(max(rowSums(otu)), big.mark = ",")),
  "- relative elevation derived from G1 reference plots",
  "",
  "Habitat comparison:",
  sprintf("- rarefied Shannon LMM p = %.3g; 0.99-coverage sqrt(q1) LMM p = %.3g",
          habitat_p, coverage_habitat_table$p),
  "",
  "Coverage:",
  sprintf("- median observed coverage: litter %.3f, soil %.3f",
          coverage_diagnostics$median_observed_coverage[coverage_diagnostics$community == "All bacteria" & coverage_diagnostics$habitat == "litter"],
          coverage_diagnostics$median_observed_coverage[coverage_diagnostics$community == "All bacteria" & coverage_diagnostics$habitat == "soil"]),
  sprintf("- sensitivity coverage reached by all libraries within 2x reads: %s",
          format(sensitivity_coverage)),
  "",
  "Primary geographic LMMs (relative elevation, additive; p and BH-adjusted p):",
  sprintf("- litter: elevation p = %.4f (BH %.4f); latitude p = %.4f (BH %.4f)",
          geo_value("litter", "relative", "elevation_p"),
          geo_value("litter", "relative", "elevation_p_BH"),
          geo_value("litter", "relative", "latitude_p"),
          geo_value("litter", "relative", "latitude_p_BH")),
  sprintf("- soil: elevation p = %.4f (BH %.4f); latitude p = %.4f (BH %.4f)",
          geo_value("soil", "relative", "elevation_p"),
          geo_value("soil", "relative", "elevation_p_BH"),
          geo_value("soil", "relative", "latitude_p"),
          geo_value("soil", "relative", "latitude_p_BH")),
  sprintf("- raw elevation: litter latitude p = %.4f; soil elevation p = %.4f",
          geo_value("litter", "raw", "latitude_p"), geo_value("soil", "raw", "elevation_p")),
  "",
  "Residual analyses (litter latitude p / soil elevation p):",
  unlist(lapply(c("relative", "raw"), function(elevation_label) {
    vapply(names(environment_sets), function(set_name) {
      sprintf("- %s elevation, %s: litter latitude p = %.3f; soil elevation p = %.3f",
              elevation_label, set_name,
              residual_value("litter", elevation_label, set_name, "residual_latitude_p"),
              residual_value("soil", elevation_label, set_name, "residual_elevation_p"))
    }, character(1))
  })),
  "",
  "Random Forest (rarefied Shannon diversity):",
  vapply(seq_len(nrow(rf_performance)), function(i) {
    sprintf("- %s: OOB R2 = %.3f; canopy axes in top 10 = %d",
            rf_performance$habitat[i], rf_performance$OOB_R2[i],
            rf_performance$top10_canopy[i])
  }, character(1)),
  "",
  sprintf("Canopy PCA, first eight axes: %.1f%% of the variance",
          100 * sum(canopy_variance)),
  sprintf("Ca-pH Pearson correlation (one row per plot): r = %.3f", ca_ph_correlation)
)
writeLines(report_lines, file.path(output_dir, "RUN_SUMMARY.txt"), useBytes = TRUE)

message("Analysis complete.")
message("See: ", file.path(output_dir, "RUN_SUMMARY.txt"))
