#!/usr/bin/env Rscript

# Reproducible analysis for:
# "Bacterial Diversity Patterns and Their Environmental Predictors in
#  European Mountain Forest Habitats"
#
# Starting point
# --------------
# This script reproduces the statistical analyses and figures using the
# processed OTU table and associated environmental metadata in `data/`.
# Raw sequence processing was performed in SEED as described in the manuscript
# and is therefore outside the scope of this R analysis script.
#
# Run from the repository root:
#   Rscript reproduce_analysis.R
#
# Optional arguments:
#   --output=/path/to/output_directory
#   --skip-pdp       Skip the partial-dependence panels (faster test run)
#   --corrected-pca  Include Betula pendula in the canopy PCA sensitivity run
#
# Required R packages:
#   iNEXT, lme4, lmerTest, MuMIn, randomForest, vegan, ggplot2,
#   patchwork, pdp, corrplot

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

# Use the project-local package library when present, without requiring it.
project_library <- file.path(script_dir, ".Rlib")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

output_arg <- grep("^--output=", command_args, value = TRUE)
output_dir <- if (length(output_arg) == 1L) {
  normalizePath(sub("^--output=", "", output_arg), mustWork = FALSE)
} else {
  file.path(script_dir, "results")
}

run_pdp <- !"--skip-pdp" %in% command_args
reproduce_archived_pca <- !"--corrected-pca" %in% command_args

data_dir <- file.path(script_dir, "data")
phyloseq_file <- file.path(data_dir, "phyloseq_16S_OTU.rds")
canopy_file <- file.path(data_dir, "canopy_vegetation.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
model_dir <- file.path(output_dir, "models")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "iNEXT", "lme4", "lmerTest", "MuMIn", "randomForest", "vegan",
  "ggplot2", "patchwork", "pdp", "corrplot"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "), "\n",
    "Install them with:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(iNEXT)
  library(lme4)
  library(lmerTest)
  library(MuMIn)
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

safe_r2 <- function(model) {
  result <- suppressWarnings(MuMIn::r.squaredGLMM(model))
  c(marginal = unname(result[1, 1]), conditional = unname(result[1, 2]))
}

extract_term_p <- function(model, variable) {
  coefficients <- coef(summary(model))
  p_column <- grep("^Pr\\(", colnames(coefficients), value = TRUE)
  if (length(p_column) == 0L) return(NA_real_)

  direct <- variable
  polynomial <- paste0("poly(", variable, ", 2)1")
  if (direct %in% rownames(coefficients)) {
    return(unname(coefficients[direct, p_column[[1]]]))
  }
  if (polynomial %in% rownames(coefficients)) {
    return(unname(coefficients[polynomial, p_column[[1]]]))
  }
  NA_real_
}

extract_interaction_p <- function(model, elevation_variable) {
  coefficients <- coef(summary(model))
  p_column <- grep("^Pr\\(", colnames(coefficients), value = TRUE)
  if (length(p_column) == 0L) return(NA_real_)
  interaction_rows <- grepl(":", rownames(coefficients)) &
    grepl(elevation_variable, rownames(coefficients), fixed = TRUE) &
    grepl("Latitude", rownames(coefficients), fixed = TRUE)
  if (!any(interaction_rows)) return(NA_real_)
  min(coefficients[interaction_rows, p_column[[1]]], na.rm = TRUE)
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

fit_candidate_lmms <- function(data, response = "sqrt_q1",
                               elevation_variable = "Elevation3",
                               interaction = FALSE) {
  formulas <- candidate_formulas(response, elevation_variable, interaction)
  ml_models <- lapply(formulas, function(formula) {
    suppressMessages(suppressWarnings(lmer(formula, data = data, REML = FALSE)))
  })
  aic_values <- vapply(ml_models, AIC, numeric(1))
  best_name <- names(which.min(aic_values))
  # Refit directly instead of update(); direct refitting is compatible with
  # both the archived lme4 series and lme4 2.x.
  best_reml <- suppressMessages(suppressWarnings(lmer(
    formulas[[best_name]], data = data, REML = TRUE
  )))

  list(
    formulas = formulas,
    ml_models = ml_models,
    aic = aic_values,
    best_name = best_name,
    best_model = best_reml,
    elevation_p = extract_term_p(best_reml, elevation_variable),
    latitude_p = extract_term_p(best_reml, "Latitude"),
    interaction_p = if (interaction) {
      extract_interaction_p(best_reml, elevation_variable)
    } else {
      NA_real_
    },
    singular = lme4::isSingular(best_reml, tol = 1e-4),
    r2 = safe_r2(best_reml)
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
    selected = names(fit$aic) == fit$best_name
  )
}

model_summary_row <- function(fit, habitat, elevation_type, interaction,
                              n_samples) {
  data.frame(
    habitat = habitat,
    elevation_type = elevation_type,
    interaction = interaction,
    n_samples = n_samples,
    selected_model = fit$best_name,
    elevation_p = fit$elevation_p,
    latitude_p = fit$latitude_p,
    interaction_p = fit$interaction_p,
    singular_fit = fit$singular,
    R2_marginal = fit$r2[["marginal"]],
    R2_conditional = fit$r2[["conditional"]]
  )
}

estimate_coverage_q1 <- function(sample_by_taxon, metadata,
                                 coverage = 0.99,
                                 prune_taxa_total = 1,
                                 prune_sample_total = 1) {
  taxon_by_sample <- t(sample_by_taxon)
  taxon_by_sample <- taxon_by_sample[
    rowSums(taxon_by_sample) > prune_taxa_total, , drop = FALSE
  ]
  taxon_by_sample <- taxon_by_sample[
    , colSums(taxon_by_sample) > prune_sample_total, drop = FALSE
  ]
  # iNEXT <=2 returned one column per Hill order, whereas iNEXT >=3 returns a
  # long table. nboot=0 avoids unnecessary confidence-interval bootstrapping
  # and is the current-version equivalent of conf=NULL in the archived code.
  if (utils::packageVersion("iNEXT") >= "3.0.0") {
    estimates <- suppressWarnings(iNEXT::estimateD(
      taxon_by_sample, datatype = "abundance", base = "coverage",
      level = coverage, nboot = 0
    ))
    q1_rows <- estimates$Order.q == 1 & estimates$SC >= coverage - 1e-8
    q1 <- estimates$qD[q1_rows]
    names(q1) <- estimates$Assemblage[q1_rows]
  } else {
    estimates <- suppressWarnings(iNEXT::estimateD(
      taxon_by_sample, datatype = "abundance", base = "coverage",
      level = coverage, conf = NULL
    ))
    q1 <- estimates[["q = 1"]]
    names(q1) <- colnames(taxon_by_sample)
  }
  q1 <- q1[!is.na(q1)]
  sample_ids <- intersect(names(q1), rownames(metadata))
  result <- metadata[sample_ids, , drop = FALSE]
  result$sample_id <- rownames(result)
  result$q1 <- unname(q1[sample_ids])
  result$sqrt_q1 <- sqrt(result$q1)
  result
}

format_p <- function(x) {
  ifelse(is.na(x), "p = NA",
         ifelse(x < 0.001, "p < 0.001", sprintf("p = %.3f", x)))
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
metadata$phylum_sample_depth <- rowSums(otu)

dataset_summary <- data.frame(
  statistic = c(
    "samples", "litter_samples", "soil_samples", "taxa", "total_reads",
    "minimum_reads", "maximum_reads", "mean_reads"
  ),
  value = c(
    nrow(otu), sum(metadata$habitat == "litter"),
    sum(metadata$habitat == "soil"), ncol(otu), sum(otu),
    min(rowSums(otu)), max(rowSums(otu)), mean(rowSums(otu))
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
# 3. Figure 1: rarefied Shannon diversity by habitat
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
habitat_p <- habitat_anova[[grep("Pr\\(", names(habitat_anova), value = TRUE)[1]]][1]

habitat_statistics <- do.call(rbind, lapply(split(habitat_data, habitat_data$habitat),
                                            function(x) {
  data.frame(
    habitat = as.character(x$habitat[1]), n = nrow(x),
    mean = mean(x$shannon), sd = sd(x$shannon),
    median = median(x$shannon), minimum = min(x$shannon), maximum = max(x$shannon)
  )
}))
habitat_statistics$LMM_p <- habitat_p
write_table(habitat_statistics, "02_habitat_shannon_statistics.csv")

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

ggsave(file.path(figure_dir, "Fig1_reproduced.pdf"), fig1,
       width = 6, height = 5, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig1_reproduced.png"), fig1,
       width = 6, height = 5, dpi = 220)

# -----------------------------------------------------------------------------
# 4. Canopy vegetation PCA (used by the Random Forest analysis)
# -----------------------------------------------------------------------------

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
canopy_species <- canopy_species[!grepl("1", rownames(canopy_species)), , drop = FALSE]
canopy_species <- canopy_species[, colSums(canopy_species, na.rm = TRUE) > 0,
                                 drop = FALSE]

# IMPORTANT REPRODUCIBILITY NOTE
# The archived PCA scripts used `species_data[, -1]` after ID_sample had already
# been removed. Consequently, the first species column (Betula pendula) was not
# included in the published PCA. TRUE reproduces the reported 52.9% variance;
# use --corrected-pca to include that species as a transparent sensitivity run.
archived_excluded_species <- NA_character_
if (reproduce_archived_pca) {
  archived_excluded_species <- names(canopy_species)[1]
  canopy_matrix <- as.matrix(canopy_species[, -1, drop = FALSE])
} else {
  canopy_matrix <- as.matrix(canopy_species)
}

canopy_pca <- vegan::rda(canopy_matrix, scale = TRUE)
canopy_scores <- as.data.frame(scores(canopy_pca, display = "sites", choices = 1:8))
canopy_loadings <- as.data.frame(scores(canopy_pca, display = "species", choices = 1:8))
canopy_variance <- summary(canopy_pca)$cont$importance[
  "Proportion Explained", 1:8
]

canopy_variance_table <- data.frame(
  axis = paste0("PC", 1:8),
  proportion_explained = as.numeric(canopy_variance),
  cumulative_first_eight = sum(canopy_variance),
  archived_reproduction = reproduce_archived_pca,
  archived_excluded_species = archived_excluded_species
)
write_table(canopy_variance_table, "03_canopy_pca_variance.csv")

canopy_scores_output <- data.frame(plot = rownames(canopy_scores), canopy_scores,
                                   row.names = NULL)
write_table(canopy_scores_output, "04_canopy_pca_scores.csv")

canopy_loadings_output <- data.frame(species = rownames(canopy_loadings),
                                     canopy_loadings, row.names = NULL)
write_table(canopy_loadings_output, "05_canopy_pca_loadings.csv")

# Compare recomputed axes with the archived axes in the final RDS. PCA signs are
# arbitrary, so absolute correlations are the meaningful diagnostic.
unique_metadata <- metadata[!duplicated(metadata$Plot), , drop = FALSE]
rownames(unique_metadata) <- unique_metadata$Plot
common_pca_plots <- intersect(rownames(canopy_scores), rownames(unique_metadata))
archived_axis_names <- paste0("PC", 1:8, ".x")
axis_correlations <- cor(
  canopy_scores[common_pca_plots, , drop = FALSE],
  unique_metadata[common_pca_plots, archived_axis_names, drop = FALSE],
  use = "pairwise.complete.obs"
)
axis_correlation_table <- data.frame(
  recomputed_axis = rownames(axis_correlations),
  archived_axis = colnames(axis_correlations)[max.col(abs(axis_correlations), ties.method = "first")],
  absolute_correlation = apply(abs(axis_correlations), 1, max)
)
write_table(axis_correlation_table, "06_canopy_pca_archived_axis_correlations.csv")

# -----------------------------------------------------------------------------
# 5. Figure 2: geographic alpha-diversity GLMMs and raw-elevation sensitivity
# -----------------------------------------------------------------------------

coverage_data <- list()
for (habitat_name in c("litter", "soil")) {
  sample_ids <- rownames(metadata)[metadata$habitat == habitat_name]
  coverage_data[[habitat_name]] <- estimate_coverage_q1(
    otu[sample_ids, , drop = FALSE], metadata, coverage = 0.99,
    prune_taxa_total = 1, prune_sample_total = 1
  )
  coverage_data[[habitat_name]]$habitat <- habitat_name
}

coverage_all <- do.call(rbind, coverage_data)
rownames(coverage_all) <- coverage_all$sample_id
write_table(coverage_all[, c("sample_id", "habitat", "q1", "sqrt_q1",
                            "Elevation3", "Elevation", "Latitude", "Location",
                            "PermGroup", "Mountain_range")],
            "07_coverage_q1_overall.csv")

geographic_fits <- list()
geographic_model_rows <- list()
geographic_aic_rows <- list()
for (habitat_name in names(coverage_data)) {
  dat <- coverage_data[[habitat_name]]
  elevation_variables <- c(relative = "Elevation3", raw = "Elevation")
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
write_table(geographic_model_table, "08_geographic_GLMM_results.csv")
write_table(geographic_aic_table, "09_geographic_GLMM_AIC_selection.csv")

relative_additive_fits <- list(
  litter = geographic_fits[["litter_relative_FALSE"]],
  soil = geographic_fits[["soil_relative_FALSE"]]
)

make_geographic_panel <- function(data, habitat_name, gradient,
                                  fit, show_legend = FALSE) {
  p_value <- if (gradient == "Elevation3") fit$elevation_p else fit$latitude_p
  label <- NULL
  if (!is.na(p_value) && p_value < 0.05) {
    label <- sprintf(
      "R²m = %.3f\nR²c = %.3f\n%s",
      fit$r2[["marginal"]], fit$r2[["conditional"]], format_p(p_value)
    )
  }

  panel <- ggplot(data, aes(x = .data[[gradient]], y = q1,
                            color = Mountain_range)) +
    geom_point(alpha = 0.72, size = 1.8) +
    scale_color_brewer(palette = "Paired") +
    labs(
      x = if (gradient == "Elevation3") "Elevation" else "Latitude",
      y = "Hill number q=1",
      color = "Mountain Range"
    ) +
    theme_figures(10) +
    theme(legend.position = if (show_legend) "right" else "none")

  if (!is.na(p_value) && p_value < 0.05) {
    polynomial <- if (gradient == "Elevation3") {
      fit$best_name %in% c("M1", "M2")
    } else {
      fit$best_name %in% c("M1", "M3")
    }
    panel <- panel + geom_smooth(
      method = "lm", formula = if (polynomial) y ~ poly(x, 2) else y ~ x,
      color = "black", fill = "grey65", linewidth = 0.7
    ) + annotate(
      "text", x = if (gradient == "Elevation3") min(data[[gradient]]) else max(data[[gradient]]),
      y = max(data$q1), label = label,
      hjust = if (gradient == "Elevation3") 0 else 1, vjust = 1, size = 3.5
    )
  }
  panel
}

fig2a <- make_geographic_panel(coverage_data$litter, "litter", "Elevation3",
                               relative_additive_fits$litter)
fig2b <- make_geographic_panel(coverage_data$litter, "litter", "Latitude",
                               relative_additive_fits$litter, show_legend = TRUE)
fig2c <- make_geographic_panel(coverage_data$soil, "soil", "Elevation3",
                               relative_additive_fits$soil)
fig2d <- make_geographic_panel(coverage_data$soil, "soil", "Latitude",
                               relative_additive_fits$soil)

fig2 <- ((fig2a + ggtitle("a) Elevation - Litter")) |
           (fig2b + ggtitle("b) Latitude - Litter"))) /
  ((fig2c + ggtitle("c) Elevation - Soil")) |
     (fig2d + ggtitle("d) Latitude - Soil"))) +
  plot_annotation(theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(figure_dir, "Fig2_reproduced.pdf"), fig2,
       width = 13.5, height = 9.5, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig2_reproduced.png"), fig2,
       width = 13.5, height = 9.5, dpi = 180)

# -----------------------------------------------------------------------------
# 6. Residual analyses: climate, soil, and combined environmental models
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

run_residual_analysis <- function(data, response = "sqrt_q1") {
  rows <- list()
  models <- list()
  for (set_name in names(environment_sets)) {
    variables <- environment_sets[[set_name]]
    complete <- complete.cases(data[, c(response, variables, "Elevation3",
                                         "Latitude", "PermGroup", "Location")])
    model_data <- data[complete, , drop = FALSE]
    environment_formula <- as.formula(sprintf(
      "%s ~ %s", response, paste(variables, collapse = " + ")
    ))
    environment_model <- lm(environment_formula, data = model_data)
    model_data$environment_residual <- residuals(environment_model)
    geographic_fit <- fit_candidate_lmms(
      model_data, response = "environment_residual",
      elevation_variable = "Elevation3", interaction = FALSE
    )
    rows[[set_name]] <- data.frame(
      environmental_set = set_name,
      n_samples = nrow(model_data),
      environmental_R2 = summary(environment_model)$r.squared,
      selected_geographic_model = geographic_fit$best_name,
      singular_fit = geographic_fit$singular,
      residual_elevation_p = geographic_fit$elevation_p,
      residual_latitude_p = geographic_fit$latitude_p
    )
    models[[set_name]] <- list(
      environment_model = environment_model,
      geographic_fit = geographic_fit,
      data = model_data
    )
  }
  list(table = do.call(rbind, rows), models = models)
}

overall_residual_results <- list()
overall_residual_rows <- list()
for (habitat_name in names(coverage_data)) {
  result <- run_residual_analysis(coverage_data[[habitat_name]])
  result$table$habitat <- habitat_name
  overall_residual_results[[habitat_name]] <- result
  overall_residual_rows[[habitat_name]] <- result$table
}
overall_residual_table <- do.call(rbind, overall_residual_rows)
overall_residual_table <- overall_residual_table[
  , c("habitat", setdiff(names(overall_residual_table), "habitat"))
]
write_table(overall_residual_table, "10_overall_residual_analysis.csv")

# Environmental correlation matrix for Fig. S1. The figure shows the complete
# candidate set before collinearity filtering, ordered as climate, canopy, and
# soil. Environmental measurements occur twice in sample_data (paired litter
# and soil); use one row per plot here.
rf_climate_candidates <- c(
  "Mean_Diurnal_Range", "MAT", "MAP", "Temperature_Seasonality",
  "Max_Temp_Warmest_Month", "Min_Temp_Coldest_Month",
  "Temperature_Annual_Range", "Mean_Temp_Wettest_Quarter",
  "Mean_Temp_Driest_Quarter", "Mean_Temp_Warmest_Quarter",
  "Mean_Temp_Coldest_Quarter", "Precipitation_Wettest_Month",
  "Precipitation_Driest_Month", "Precipitation_Seasonality",
  "Precipitation_Wettest_Quarter", "Precipitation_Driest_Quarter",
  "Precipitation_Warmest_Quarter", "Precipitation_Coldest_Quarter"
)
rf_soil_candidates <- c("c", "n", "Ca", "pH", "p_MEL")
rf_canopy_candidates <- archived_axis_names
correlation_variables <- c(
  rf_climate_candidates, rf_canopy_candidates, rf_soil_candidates
)
correlation_metadata <- metadata[!duplicated(metadata$Plot), , drop = FALSE]
correlation_data <- correlation_metadata[, correlation_variables, drop = FALSE]
correlation_matrix <- cor(correlation_data, use = "pairwise.complete.obs")
write_table(data.frame(variable = rownames(correlation_matrix), correlation_matrix,
                       check.names = FALSE),
            "11_environmental_correlation_matrix.csv")

# Also retain the matrix for the reduced predictor set used in the final
# residual and Random Forest analyses.
selected_correlation_variables <- c(
  climate_variables, archived_axis_names, soil_variables
)
selected_correlation_matrix <- cor(
  correlation_metadata[, selected_correlation_variables, drop = FALSE],
  use = "pairwise.complete.obs"
)
write_table(data.frame(
  variable = rownames(selected_correlation_matrix), selected_correlation_matrix,
  check.names = FALSE
), "11b_selected_predictor_correlation_matrix.csv")

# Use the same corrplot call as the archived OTU.R analysis.
grDevices::pdf(file.path(figure_dir, "FigS1_correlation_reproduced.pdf"),
               width = 11, height = 9)
corrplot::corrplot(correlation_matrix, method = "circle", type = "upper",
                   tl.col = "black", tl.cex = 0.7)
grDevices::dev.off()

grDevices::png(file.path(figure_dir, "FigS1_correlation_reproduced.png"),
               width = 11, height = 9, units = "in", res = 220)
corrplot::corrplot(correlation_matrix, method = "circle", type = "upper",
                   tl.col = "black", tl.cex = 0.7)
grDevices::dev.off()

# -----------------------------------------------------------------------------
# 7. Phylum-level coverage, model selection, residual analyses, Figures 3-4
# -----------------------------------------------------------------------------

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
), "12_phylum_mean_relative_abundance.csv")

exclude_litter_phyla <- c("Chloroflexota", "Gemmatimonadota", "Patescibacteriota")
phylum_analysis <- list()
phylum_rows <- list()

for (habitat_name in c("soil", "litter")) {
  if (habitat_name == "litter") {
    phyla_to_run <- setdiff(major_phyla, exclude_litter_phyla)
  } else {
    phyla_to_run <- major_phyla
  }
  habitat_samples <- rownames(metadata)[metadata$habitat == habitat_name]
  for (phylum_name in phyla_to_run) {
    taxon_keep <- phylum == phylum_name
    phylum_data <- estimate_coverage_q1(
      otu[habitat_samples, taxon_keep, drop = FALSE], metadata,
      coverage = 0.99, prune_taxa_total = 1, prune_sample_total = 1
    )
    phylum_data$habitat <- habitat_name
    fit <- fit_candidate_lmms(phylum_data, response = "sqrt_q1",
                              elevation_variable = "Elevation3",
                              interaction = FALSE)
    key <- paste(phylum_name, habitat_name, sep = "_")
    phylum_analysis[[key]] <- list(data = phylum_data, fit = fit)
    phylum_rows[[key]] <- data.frame(
      phylum = phylum_name,
      habitat = habitat_name,
      n_samples = nrow(phylum_data),
      coverage = 0.99,
      selected_model = fit$best_name,
      AIC = min(fit$aic),
      singular_fit = fit$singular,
      elevation_p = fit$elevation_p,
      latitude_p = fit$latitude_p,
      R2_marginal = fit$r2[["marginal"]],
      R2_conditional = fit$r2[["conditional"]]
    )
  }
}

phylum_results <- do.call(rbind, phylum_rows)
phylum_results$elevation_p_BH <- p.adjust(phylum_results$elevation_p, method = "BH")
phylum_results$latitude_p_BH <- p.adjust(phylum_results$latitude_p, method = "BH")
write_table(phylum_results, "13_phylum_GLMM_results.csv")

# Residual analyses are required only for phylum-habitat combinations with at
# least one significant geographic relationship in the main GLMM analysis.
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
  phylum_residual_rows[[key]] <- residual_result$table
}
phylum_residual_table <- do.call(rbind, phylum_residual_rows)
phylum_residual_table <- phylum_residual_table[
  , c("phylum", "habitat", setdiff(names(phylum_residual_table),
                                    c("phylum", "habitat")))
]
write_table(phylum_residual_table, "14_phylum_residual_analysis.csv")

make_phylum_figure <- function(gradient = c("Elevation3", "Latitude")) {
  gradient <- match.arg(gradient)
  p_column <- if (gradient == "Elevation3") "elevation_p_BH" else "latitude_p_BH"
  significant_phyla <- unique(phylum_results$phylum[
    !is.na(phylum_results[[p_column]]) & phylum_results[[p_column]] < 0.05
  ])
  plot_rows <- list()
  annotation_rows <- list()

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
    dat$significant <- !is.na(result_row[[p_column]]) && result_row[[p_column]] < 0.05
    dat$model <- result_row$selected_model
    plot_rows[[key]] <- dat

    if (dat$significant[1]) {
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
          format_p(result_row[[p_column]])
        )
      )
    }
  }

  plot_data <- do.call(rbind, plot_rows)
  annotations <- do.call(rbind, annotation_rows)
  plot_data$habitat <- factor(plot_data$habitat, levels = c("litter", "soil"))

  figure <- ggplot(plot_data, aes(x = .data[[gradient]], y = q1,
                                  color = habitat)) +
    geom_point(aes(alpha = significant), size = 1.5) +
    scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.25), guide = "none") +
    scale_color_manual(values = c(litter = "#01665e", soil = "#8c510a"),
                       labels = c(litter = "Litter", soil = "Soil")) +
    facet_wrap(~phylum, scales = "free_y", ncol = 3) +
    labs(
      x = if (gradient == "Elevation3") "Elevation" else "Latitude (°)",
      y = "Hill number q=1", color = "Sample Type"
    ) +
    theme_figures(10) +
    theme(legend.position = "bottom")

  for (key in names(phylum_analysis)) {
    entry <- phylum_analysis[[key]]
    phylum_name <- sub("_(soil|litter)$", "", key)
    habitat_name <- sub("^.*_", "", key)
    if (!phylum_name %in% significant_phyla) next
    result_row <- phylum_results[
      phylum_results$phylum == phylum_name & phylum_results$habitat == habitat_name,
    ]
    if (is.na(result_row[[p_column]]) || result_row[[p_column]] >= 0.05) next
    dat <- entry$data
    dat$phylum <- phylum_name
    dat$habitat <- habitat_name
    polynomial <- if (gradient == "Elevation3") {
      result_row$selected_model %in% c("M1", "M2")
    } else {
      result_row$selected_model %in% c("M1", "M3")
    }
    figure <- figure + geom_smooth(
      data = dat, aes(x = .data[[gradient]], y = q1, color = habitat),
      method = "lm", formula = if (polynomial) y ~ poly(x, 2) else y ~ x,
      se = TRUE, linewidth = 0.8, inherit.aes = FALSE, show.legend = FALSE
    )
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

fig3 <- make_phylum_figure("Elevation3")
fig4 <- make_phylum_figure("Latitude")
ggsave(file.path(figure_dir, "Fig3_reproduced.pdf"), fig3,
       width = 13, height = 9, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig3_reproduced.png"), fig3,
       width = 13, height = 9, dpi = 180)
ggsave(file.path(figure_dir, "Fig4_reproduced.pdf"), fig4,
       width = 13, height = 9, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig4_reproduced.png"), fig4,
       width = 13, height = 9, dpi = 180)

# -----------------------------------------------------------------------------
# 8. Random Forest tuning, OOB performance, importance, PDPs, Figure 5
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
rf_data$shannon <- shannon_index(otu)[rownames(rf_data)]

rf_climate_selected <- reduce_collinearity(rf_data, rf_climate_candidates)
rf_soil_selected <- reduce_collinearity(rf_data, rf_soil_candidates)
rf_canopy_selected <- reduce_collinearity(rf_data, rf_canopy_candidates)
rf_predictors <- c(rf_climate_selected, rf_soil_selected, rf_canopy_selected)

rf_predictor_table <- data.frame(
  variable = rf_predictors,
  category = ifelse(
    rf_predictors %in% rf_climate_selected, "Climate",
    ifelse(rf_predictors %in% rf_soil_selected, "Soil", "Canopy")
  )
)
write_table(rf_predictor_table, "15_random_forest_selected_predictors.csv")

run_random_forest <- function(data, habitat_name) {
  model_data <- data[complete.cases(data[, c("shannon", rf_predictors)]),
                     c("shannon", rf_predictors), drop = FALSE]
  mtry_values <- seq(
    floor(sqrt(length(rf_predictors))),
    min(ceiling(length(rf_predictors) / 3), length(rf_predictors)), by = 1
  )
  tuning_rows <- list()
  for (mtry_value in mtry_values) {
    set.seed(123)
    tuning_model <- randomForest(
      shannon ~ ., data = model_data, mtry = mtry_value,
      ntree = 500, importance = FALSE
    )
    tuning_rows[[as.character(mtry_value)]] <- data.frame(
      habitat = habitat_name, mtry = mtry_value,
      OOB_MSE = tail(tuning_model$mse, 1)
    )
  }
  tuning <- do.call(rbind, tuning_rows)
  best_mtry <- tuning$mtry[which.min(tuning$OOB_MSE)]

  set.seed(456)
  model <- randomForest(
    shannon ~ ., data = model_data, mtry = best_mtry,
    ntree = 1000, importance = TRUE
  )
  importance_matrix <- randomForest::importance(model, type = 1)
  importance <- data.frame(
    variable = rownames(importance_matrix),
    importance = importance_matrix[, 1],
    category = ifelse(
      rownames(importance_matrix) %in% rf_climate_selected, "Climate",
      ifelse(rownames(importance_matrix) %in% rf_soil_selected, "Soil", "Canopy")
    ),
    habitat = habitat_name
  )
  importance <- importance[order(importance$importance, decreasing = TRUE), ]

  oob_mse <- tail(model$mse, 1)
  performance <- data.frame(
    habitat = habitat_name,
    n_samples = nrow(model_data),
    n_predictors = length(rf_predictors),
    ntree = model$ntree,
    mtry = model$mtry,
    OOB_MSE = oob_mse,
    OOB_RMSE = sqrt(oob_mse),
    OOB_R2 = 1 - oob_mse / var(model$y),
    percent_variance_explained = 100 * (1 - oob_mse / var(model$y))
  )

  list(model = model, model_data = model_data, tuning = tuning,
       importance = importance, performance = performance)
}

rf_results <- list()
for (habitat_name in c("litter", "soil")) {
  rf_results[[habitat_name]] <- run_random_forest(
    rf_data[rf_data$habitat == habitat_name, , drop = FALSE], habitat_name
  )
}

rf_tuning <- do.call(rbind, lapply(rf_results, `[[`, "tuning"))
rf_importance <- do.call(rbind, lapply(rf_results, `[[`, "importance"))
rf_performance <- do.call(rbind, lapply(rf_results, `[[`, "performance"))
write_table(rf_tuning, "16_random_forest_tuning.csv")
write_table(rf_performance, "17_random_forest_performance.csv")
write_table(rf_importance, "18_random_forest_importance.csv")

variable_labels <- c(
  PC1.x = "PC1 canopy", PC2.x = "PC2 canopy", PC3.x = "PC3 canopy",
  PC4.x = "PC4 canopy", PC5.x = "PC5 canopy", PC6.x = "PC6 canopy",
  PC7.x = "PC7 canopy", PC8.x = "PC8 canopy",
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

ggsave(file.path(figure_dir, "Fig5_reproduced.pdf"), fig5,
       width = 11, height = 8, device = grDevices::pdf)
ggsave(file.path(figure_dir, "Fig5_reproduced.png"), fig5,
       width = 11, height = 8, dpi = 200)

if (run_pdp) {
  make_pdp_plot <- function(result, variable, habitat_name) {
    partial_result <- pdp::partial(
      result$model, pred.var = variable, train = result$model_data,
      progress = "none"
    )
    autoplot(partial_result, rug = TRUE, train = result$model_data) +
      labs(x = label_variable(variable), y = "Partial effect") +
      theme_minimal(base_size = 9)
  }

  for (habitat_name in c("litter", "soil")) {
    variables <- head(rf_results[[habitat_name]]$importance$variable, 10)
    plots <- lapply(variables, function(variable) {
      make_pdp_plot(rf_results[[habitat_name]], variable, habitat_name)
    })
    pdp_figure <- wrap_plots(plots, ncol = 2) +
      plot_annotation(title = paste("Partial dependence -", tools::toTitleCase(habitat_name)))
    figure_number <- if (habitat_name == "litter") "S2" else "S3"
    ggsave(file.path(figure_dir,
                     paste0("Fig", figure_number, "_PDP_", habitat_name,
                            "_reproduced.pdf")),
           pdp_figure, width = 8, height = 10, device = grDevices::pdf)
  }
}

# Species loadings for canopy axes appearing among either habitat's top ten.
important_canopy_axes <- unique(sub("\\.x$", "", rf_top10$variable[
  rf_top10$category == "Canopy"
]))
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
write_table(important_loadings, "19_canopy_top_species_loadings.csv")

loading_plots <- lapply(split(important_loadings, important_loadings$axis), function(x) {
  x$species <- factor(x$species, levels = x$species[order(x$absolute_loading)])
  axis_name <- unique(x$axis)
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
ggsave(file.path(figure_dir, "FigS4_canopy_loadings_reproduced.pdf"), fig_s4,
       width = 10, height = max(6, ceiling(length(loading_plots) / 2) * 3.5),
       device = grDevices::pdf)
ggsave(file.path(figure_dir, "FigS4_canopy_loadings_reproduced.png"), fig_s4,
       width = 10, height = max(6, ceiling(length(loading_plots) / 2) * 3.5),
       dpi = 220)

# -----------------------------------------------------------------------------
# 9. Save fitted models and run summary
# -----------------------------------------------------------------------------

get_geo <- function(habitat_name, field) {
  geographic_model_table[
    geographic_model_table$habitat == habitat_name &
      geographic_model_table$elevation_type == "relative" &
      !geographic_model_table$interaction,
    field
  ][1]
}

saveRDS(
  list(
    habitat_model = habitat_model,
    geographic_fits = geographic_fits,
    residual_results = overall_residual_results,
    phylum_analysis = phylum_analysis,
    rf_results = rf_results,
    canopy_pca = canopy_pca
  ),
  file.path(model_dir, "fitted_models.rds")
)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

report_lines <- c(
  "REPRODUCIBILITY RUN SUMMARY",
  "===========================",
  paste("Run date:", as.character(Sys.time())),
  paste("Input:", phyloseq_file),
  paste("Output:", output_dir),
  "",
  "Dataset:",
  sprintf("- %d samples (%d litter, %d soil)", nrow(otu),
          sum(metadata$habitat == "litter"), sum(metadata$habitat == "soil")),
  sprintf("- %s retained reads; depth %s-%s", format(sum(otu), big.mark = ","),
          format(min(rowSums(otu)), big.mark = ","),
          format(max(rowSums(otu)), big.mark = ",")),
  "",
  "Main geographic GLMMs:",
  sprintf("- litter, relative elevation: elevation p = %.6f; latitude p = %.6f",
          get_geo("litter", "elevation_p"), get_geo("litter", "latitude_p")),
  sprintf("- soil, relative elevation: elevation p = %.6f; latitude p = %.6f",
          get_geo("soil", "elevation_p"), get_geo("soil", "latitude_p")),
  sprintf("- litter, raw elevation sensitivity: elevation p = %.6f; latitude p = %.6f",
          geographic_model_table$elevation_p[
            geographic_model_table$habitat == "litter" &
              geographic_model_table$elevation_type == "raw" &
              !geographic_model_table$interaction][1],
          geographic_model_table$latitude_p[
            geographic_model_table$habitat == "litter" &
              geographic_model_table$elevation_type == "raw" &
              !geographic_model_table$interaction][1]),
  sprintf("- soil, raw elevation sensitivity: elevation p = %.6f; latitude p = %.6f",
          geographic_model_table$elevation_p[
            geographic_model_table$habitat == "soil" &
              geographic_model_table$elevation_type == "raw" &
              !geographic_model_table$interaction][1],
          geographic_model_table$latitude_p[
            geographic_model_table$habitat == "soil" &
              geographic_model_table$elevation_type == "raw" &
              !geographic_model_table$interaction][1]),
  "",
  "Random Forest OOB performance:",
  vapply(seq_len(nrow(rf_performance)), function(i) {
    sprintf("- %s: OOB R2 = %.3f, OOB RMSE = %.3f, mtry = %d, n = %d",
            rf_performance$habitat[i], rf_performance$OOB_R2[i],
            rf_performance$OOB_RMSE[i], rf_performance$mtry[i],
            rf_performance$n_samples[i])
  }, character(1)),
  "",
  paste0("Canopy PCA first eight axes: ", sprintf("%.1f%%", 100 * sum(canopy_variance)),
         if (reproduce_archived_pca) {
           paste0(" (archived analysis reproduced; excluded ", archived_excluded_species, ")")
         } else {
           " (corrected sensitivity analysis)"
         })
)
writeLines(report_lines, file.path(output_dir, "RUN_SUMMARY.txt"), useBytes = TRUE)

message("Analysis complete.")
message("See: ", file.path(output_dir, "RUN_SUMMARY.txt"))
