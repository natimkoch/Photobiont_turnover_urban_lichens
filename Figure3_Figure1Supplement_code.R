# =============================================================================
# Urban Photobionts — Figure 3 & Supplementary Figures (Heavy Metals)
# Manuscript: Urban photobionts manuscript (Koch et al.)
# Description: Produces scatter plots with otu-specific regression lines for
#              NDVI and heavy metals (Pb, Zn, Fe) against four ecophysiological
#              response variables. Also fits and exports linear models testing
#              otu × predictor interactions.
#
# Notes on smooth line decisions:
#   - NDVI plots: I05 smooth line excluded (non-significant / low n)
#   - Zn ~ FvFm:  A46 smooth line excluded (non-significant / low n)
#   - Heavy metal plots (Pb, Zn, Fe) for Gross/Net Assimilation and
#     Respiration: no smooth lines shown (patterns unclear across otus)
# =============================================================================

# -----------------------------------------------------------------------------
# 1. LIBRARIES
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(cowplot)
library(broom)

# -----------------------------------------------------------------------------
# 2. DATA IMPORT
# -----------------------------------------------------------------------------
# Load the merged photobiont + ecophysiology dataset.
# Assumes the file is in the current working directory; adjust path if needed.

photo_all <- read.csv(
  "photo_all_nowtp.csv",
  stringsAsFactors = TRUE,
  header           = TRUE,
  row.names        = 1,
  na.strings       = "NA",
  check.names      = FALSE
)

# -----------------------------------------------------------------------------
# 3. DOMINANT OTU ASSIGNMENT
# -----------------------------------------------------------------------------
# A sample is assigned to a otu only if that otu exceeds 70% of reads.
# Samples below this threshold are treated as NA (mixed/unassigned).

# -----------------------------------------------------------------------------
# 3. SHARED PLOT AESTHETICS
# -----------------------------------------------------------------------------

# Color palette: A46 = teal, I05 = yellow-green, third otu = gold
otu_colors <- c("#00A6AA", "#CDD986", "#C39D3B")

scale_col  <- scale_color_manual(values = otu_colors, na.value = "grey70")
scale_fill <- scale_fill_manual(values = otu_colors)

# Helper: otu-specific LM smooth, with option to exclude one otu
smooth_lm <- function(data, exclude_otu = NULL) {
  d <- data %>% filter(!is.na(max_otu))
  if (!is.null(exclude_otu)) d <- d %>% filter(max_otu != exclude_otu)
  geom_smooth(
    data      = d,
    aes(fill  = max_otu),
    method    = "lm",
    se        = TRUE,
    alpha     = 0.2,
    linewidth = 0.8
  )
}

# Helper: builds a base scatter plot (points only, no smooth)
base_plot <- function(xvar, yvar, xlabel, ylabel, ylims = NULL) {
  ggplot(photo_all,
         aes(x = .data[[xvar]], y = .data[[yvar]], color = max_otu)) +
    geom_point(size = 2.5) +
    scale_col + scale_fill +
    scale_x_continuous(name = xlabel) +
    scale_y_continuous(name = ylabel, limits = ylims) +
    theme_classic()
}

# -----------------------------------------------------------------------------
# 4. FIGURE 3 — NDVI VS. ECOPHYSIOLOGY
# -----------------------------------------------------------------------------
# I05 smooth line is excluded from all NDVI panels (see header note).

d_no_I05 <- photo_all %>% filter(!is.na(max_otu), max_otu != "I05")

ndvi_grossassim <- base_plot("ndvi10", "Assimilation_gross", "NDVI", "Gross Carbon Assimilation") +
  smooth_lm(photo_all, exclude_otu = "I05")

ndvi_assim <- base_plot("ndvi10", "Assimilation_liquid", "NDVI", "Net Carbon Assimilation") +
  smooth_lm(photo_all, exclude_otu = "I05")

# Respiration: no smooth line shown for any otu
ndvi_res <- base_plot("ndvi10", "Respiration_liquid", "NDVI", "Respiration")

ndvi_fvfm <- base_plot("ndvi10", "FvFm_liquid", "NDVI", "Fv/Fm") +
  smooth_lm(photo_all, exclude_otu = "I05")

# Figure 3: 2×2 panel
plot_grid(ndvi_grossassim, ndvi_assim, ndvi_res, ndvi_fvfm, ncol = 2, labels = "AUTO")

# -----------------------------------------------------------------------------
# 5. SUPPLEMENTARY — HEAVY METALS VS. ECOPHYSIOLOGY
# -----------------------------------------------------------------------------
# Smooth lines are shown only where indicated (see zn_fvfm below).
# All other heavy metal panels show points only.

# --- Lead (Pb) ---
pb_grossassim <- base_plot("Pb", "Assimilation_gross",  "Lead (ppm)", "Gross Carbon Assimilation")
pb_assim      <- base_plot("Pb", "Assimilation_liquid", "Lead (ppm)", "Net Carbon Assimilation",  ylims = c(0, 7.5))
pb_res        <- base_plot("Pb", "Respiration_liquid",  "Lead (ppm)", "Respiration")
pb_fvfm       <- base_plot("Pb", "FvFm_liquid",         "Lead (ppm)", "Fv/Fm", ylims = c(0, 0.9))

plot_grid(pb_grossassim, pb_assim, pb_res, pb_fvfm, ncol = 2, labels = "AUTO")

# --- Zinc (Zn) ---
zn_grossassim <- base_plot("Zn", "Assimilation_gross",  "Zinc (ppm)", "Gross Carbon Assimilation")
zn_assim      <- base_plot("Zn", "Assimilation_liquid", "Zinc (ppm)", "Net Carbon Assimilation",  ylims = c(0, 7.5))
zn_res        <- base_plot("Zn", "Respiration_liquid",  "Zinc (ppm)", "Respiration liquid")

# Zn ~ FvFm: A46 smooth excluded; I05 smooth retained
zn_fvfm <- base_plot("Zn", "FvFm_liquid", "Zinc (ppm)", "Fv/Fm", ylims = c(0, 0.9)) +
  smooth_lm(photo_all, exclude_otu = "A46")
zn_res <- base_plot("Zn", "Respiration_liquid", "Zinc (ppm)", "Respiration") +
  smooth_lm(photo_all, exclude_otu = "A46")

plot_grid(zn_grossassim, zn_assim, zn_res, zn_fvfm, ncol = 2, labels = "AUTO")

# --- Iron (Fe) ---
fe_grossassim <- base_plot("Fe", "Assimilation_gross",  "Iron (ppm)", "Gross Carbon Assimilation")
fe_assim      <- base_plot("Fe", "Assimilation_liquid", "Iron (ppm)", "Net Carbon Assimilation",  ylims = c(0, 7.5))
fe_res        <- base_plot("Fe", "Respiration_liquid",  "Iron (ppm)", "Respiration")
fe_fvfm       <- base_plot("Fe", "FvFm_liquid",         "Iron (ppm)", "Fv/Fm", ylims = c(0, 0.9))

plot_grid(fe_grossassim, fe_assim, fe_res, fe_fvfm, ncol = 2, labels = "AUTO")

# Combined 12-panel supplementary figure: rows = response variables, columns = metals
plot_grid(
  fe_grossassim, pb_grossassim, zn_grossassim,
  fe_assim,      pb_assim,      zn_assim,
  fe_res,        pb_res,        zn_res,
  fe_fvfm,       pb_fvfm,       zn_fvfm,
  ncol = 3, labels = "AUTO"
)

# -----------------------------------------------------------------------------
# 6. LINEAR MODELS — PREDICTOR × OTU INTERACTIONS
# -----------------------------------------------------------------------------
# Model structure: response ~ predictor * max_otu - 1
#   - The -1 removes the global intercept, giving each otu its own baseline.
#   - The interaction (predictor:max_otu) tests whether the slope of the
#     predictor differs significantly among otus.
#
# run_otu_models() loops over four response variables for a given predictor
# and prints each model summary to the console.

responses <- c("Assimilation_gross", "Assimilation_liquid", "Respiration_liquid", "FvFm_liquid")

run_otu_models <- function(predictor, data = photo_all) {
  lapply(responses, function(resp) {
    f <- as.formula(paste(resp, "~", predictor, "* max_otu - 1"))
    m <- lm(f, data = data)
    cat("\n===", resp, "~", predictor, "* max_otu ===\n")
    print(summary(m))
    invisible(m)
  })
}

# NDVI models (also includes null model without otu for comparison)
ndvi_models <- run_otu_models("ndvi10")

# Null NDVI models (no otu term) — used to assess overall NDVI effect
lapply(responses, function(resp) {
  f  <- as.formula(paste(resp, "~ ndvi10"))
  m  <- lm(f, data = photo_all)
  cat("\n===", resp, "~ ndvi10 (no otu) ===\n")
  print(summary(m))
})

# Heavy metal models
pb_models <- run_otu_models("Pb")
zn_models <- run_otu_models("Zn")
fe_models <- run_otu_models("Fe")

# -----------------------------------------------------------------------------
# 7. OTU-SPECIFIC SLOPES — EXPORT TO CSV
# -----------------------------------------------------------------------------
# Fits separate regressions within each dominant otu using group_modify(),
# which retains the grouping column (max_otu) in the output automatically.
# All predictor × response × otu combinations are combined into one CSV.

fit_by_otu <- function(predictor, response, data = photo_all) {
  f <- as.formula(paste(response, "~", predictor))
  data %>%
    group_by(max_otu) %>%
    group_modify(~ broom::tidy(lm(f, data = .x))) %>%
    ungroup() %>%
    mutate(predictor = predictor, response = response)
}

predictors <- c("ndvi10", "Pb", "Zn", "Fe")

all_results <- bind_rows(lapply(predictors, function(pred) {
  bind_rows(lapply(responses, function(resp) {
    fit_by_otu(pred, resp)
  }))
})) %>%
  dplyr::select(predictor, response, max_otu, term, estimate, std.error, statistic, p.value)

write.csv(all_results, "~/Documents/model_results_by_otu_all.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. SPEARMAN CORRELATION TESTS BETWEEN NDVI AND LICHEN ECOPHYSIOLOGY
# -----------------------------------------------------------------------------

cor.test(photo_all$ndvi10, photo_all$Assimilation_liquid, method = "spearman")
cor.test(photo_all$ndvi10, photo_all$Assimilation_gross, method = "spearman")
cor.test(photo_all$ndvi10, photo_all$FvFm_liquid, method = "spearman")
cor.test(photo_all$ndvi10, photo_all$Respiration_liquid, method = "spearman")
