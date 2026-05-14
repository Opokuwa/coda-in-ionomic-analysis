# =============================================================================
# Compositional Methods in Ionomic Data Analysis — PCA and Biplots
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Performs PCA on both CLR-transformed and raw-scaled elemental data.
#   Produces scores plots, loadings plots, combined biplots, and highlighted
#   aquifer plots — used in Section 3.3 (Visualization of Data Structure).
#
# INPUT:   Amphipod ionomic dataset (113 samples × 28 elements)
# OUTPUT:  Ionome_all_PCA_plots.pdf — all PCA visualizations
# =============================================================================

library(dplyr)
library(compositions)
library(ggplot2)
library(readxl)
library(grid)


# =============================================================================
# 1. DATA LOADING
# =============================================================================

# Update the path below to point to your local copy of the dataset
fn  <- "OK Springs master datasheet.xlsx"
dat <- read_excel(fn, sheet = "Pods ICP")

elem_cols <- c("Al (ug/g)", "As (ug/g)", "B (ug/g)",  "Ba (ug/g)", "Be (ug/g)",
               "Bi (ug/g)", "Ca (ug/g)", "Cd (ug/g)", "Co (ug/g)", "Cr (ug/g)",
               "Cu (ug/g)", "Fe (ug/g)", "K (ug/g)",  "Li (ug/g)", "Mg (ug/g)",
               "Mn (ug/g)", "Mo (ug/g)", "Na (ug/g)", "Ni (ug/g)", "P (ug/g)",
               "Pb (ug/g)", "S (ug/g)",  "Se (ug/g)", "Si (ug/g)", "Sr (ug/g)",
               "Tl (ug/g)", "V (ug/g)",  "Zn (ug/g)")
elem_clean <- gsub(" \\(ug/g\\)", "", elem_cols)
names(dat)[match(elem_cols, names(dat))] <- elem_clean

X_raw  <- dat[, elem_clean]
X_comp <- acomp(X_raw)
X_clr  <- clr(X_comp)


# =============================================================================
# 2. PCA
# =============================================================================
# CLR PCA: center only (CLR coordinates are already scale-invariant)
# Raw PCA: center AND scale (raw concentrations have incompatible units/ranges)

pca_clr <- prcomp(X_clr, center = TRUE, scale. = FALSE)
pca_raw <- prcomp(X_raw, center = TRUE, scale. = TRUE)

imp_clr <- summary(pca_clr)$importance
imp_raw <- summary(pca_raw)$importance

# Score data frames with grouping variables
scores_clr <- as.data.frame(pca_clr$x) %>%
  bind_cols(dat %>% dplyr::select(Aquifer, Sex, Region, Site))
scores_raw <- as.data.frame(pca_raw$x) %>%
  bind_cols(dat %>% dplyr::select(Aquifer, Sex, Region, Site))

# Loadings — retain top 50% contributors for readability
scale_factor <- 2

make_loadings <- function(pca_obj, scale_fac = 2) {
  ld <- as.data.frame(pca_obj$rotation[, 1:2])
  ld$Element <- rownames(ld)
  ld %>%
    mutate(contrib = PC1^2 + PC2^2) %>%
    filter(contrib >= quantile(contrib, 0.5)) %>%
    mutate(PC1_end = PC1 * scale_fac,
           PC2_end = PC2 * scale_fac)
}

loadings_clr_small <- make_loadings(pca_clr, scale_factor)
loadings_raw_small <- make_loadings(pca_raw, scale_factor)


# =============================================================================
# 3. GLOBAL THEME
# =============================================================================

theme_set(
  theme_bw(base_size = 18) +
    theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title    = element_text(size = 18, face = "bold"),
      axis.text     = element_text(size = 14, color = "black"),
      legend.title  = element_text(size = 14, face = "bold"),
      legend.text   = element_text(size = 12),
      legend.key.size = unit(0.7, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
)


# =============================================================================
# 4. PLOT FUNCTIONS
# =============================================================================

# --- Scores plot (colored by grouping variable) ---
plot_scores <- function(scores_df, imp, group_var, title_prefix) {
  ggplot(scores_df, aes(x = PC1, y = PC2, color = .data[[group_var]])) +
    geom_point(alpha = 0.8, size = 2.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.7) +
    labs(x     = paste0("PC1 (", round(imp[2, 1] * 100, 1), "%)"),
         y     = paste0("PC2 (", round(imp[2, 2] * 100, 1), "%)"),
         title = paste(title_prefix, "Scores by", group_var),
         color = group_var)
}

# --- Loadings plot (top-contributing elements only) ---
plot_loadings <- function(loadings_df, imp, title_prefix) {
  ggplot(loadings_df, aes(x = PC1_end, y = PC2_end, label = Element)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1_end, yend = PC2_end),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "purple", linewidth = 1, alpha = 0.9) +
    geom_text(size = 4.5, color = "black", fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    labs(x     = paste0("PC1 (", round(imp[2, 1] * 100, 1), "%)"),
         y     = paste0("PC2 (", round(imp[2, 2] * 100, 1), "%)"),
         title = paste(title_prefix, "Loadings"))
}

# --- Biplot (scores + loadings overlay) ---
make_biplot <- function(pca_obj, scores_df, group_var,
                        title_prefix = "PCA", top_prop = 0.7) {
  imp      <- summary(pca_obj)$importance
  loadings <- as.data.frame(pca_obj$rotation[, 1:2])
  loadings$Element <- rownames(loadings)
  loadings <- loadings %>%
    mutate(contrib = PC1^2 + PC2^2) %>%
    filter(contrib >= quantile(contrib, top_prop))

  # Scale loadings to fit within the scores plot window
  scale_fac <- min(
    max(abs(scores_df$PC1)) / max(abs(loadings$PC1)),
    max(abs(scores_df$PC2)) / max(abs(loadings$PC2))
  ) * 0.7
  loadings <- loadings %>%
    mutate(PC1_end = PC1 * scale_fac,
           PC2_end = PC2 * scale_fac)

  ggplot() +
    geom_point(data = scores_df,
               aes(x = PC1, y = PC2, color = .data[[group_var]]),
               alpha = 0.75, size = 2.6) +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1_end, yend = PC2_end),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "purple", linewidth = 1) +
    geom_text(data = loadings,
              aes(x = PC1_end, y = PC2_end, label = Element),
              size = 4, fontface = "bold", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    labs(x     = paste0("PC1 (", round(imp[2, 1] * 100, 1), "%)"),
         y     = paste0("PC2 (", round(imp[2, 2] * 100, 1), "%)"),
         title = paste(title_prefix, "Biplot by", group_var),
         color = group_var)
}

# --- Highlighted aquifer plot ---
# Greys out all other aquifer samples to focus on one aquifer at a time
highlight_aquifer_pca <- function(scores_df, imp, aquifer_name,
                                  title_prefix = "Ionome PCA") {
  focal  <- scores_df %>% filter(Aquifer == aquifer_name)
  others <- scores_df %>% filter(Aquifer != aquifer_name)

  ggplot() +
    geom_point(data = others, aes(PC1, PC2),
               color = "grey85", alpha = 0.4, size = 2.2) +
    geom_point(data = focal, aes(PC1, PC2, color = Site),
               size = 3, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    labs(x     = paste0("PC1 (", round(imp[2, 1] * 100, 1), "%)"),
         y     = paste0("PC2 (", round(imp[2, 2] * 100, 1), "%)"),
         title = paste0(title_prefix, ": ", aquifer_name, " highlighted"),
         color = "Site")
}


# =============================================================================
# 5. PDF OUTPUT — ALL PLOTS
# =============================================================================

group_vars <- c("Aquifer", "Sex", "Region", "Site")

pdf("Ionome_all_PCA_plots.pdf", width = 12, height = 8)

# 1. Scores plots (CLR and Raw, each grouping variable)
for (group in group_vars) {
  print(plot_scores(scores_clr, imp_clr, group, "Ionome PCA (CLR)"))
  print(plot_scores(scores_raw, imp_raw, group, "Ionome PCA (Raw)"))
}

# 2. Loadings plots
print(plot_loadings(loadings_clr_small, imp_clr, "Ionome PCA (CLR)"))
print(plot_loadings(loadings_raw_small, imp_raw, "Ionome PCA (Raw)"))

# 3. Biplots
for (group in group_vars) {
  print(make_biplot(pca_clr, scores_clr, group, "Ionome PCA (CLR)"))
  print(make_biplot(pca_raw, scores_raw, group, "Ionome PCA (Raw)"))
}

# 4. Highlighted aquifer plots
scores_clr_aq <- as.data.frame(pca_clr$x) %>%
  bind_cols(dat %>% dplyr::select(Aquifer, Site))
scores_raw_aq <- as.data.frame(pca_raw$x) %>%
  bind_cols(dat %>% dplyr::select(Aquifer, Site))

aquifers <- sort(unique(dat$Aquifer))
for (aq in aquifers) {
  print(highlight_aquifer_pca(scores_clr_aq, imp_clr, aq, "Ionome PCA (CLR)"))
  print(highlight_aquifer_pca(scores_raw_aq, imp_raw, aq, "Ionome PCA (Raw)"))
}

dev.off()
cat("Saved all PCA plots to Ionome_all_PCA_plots.pdf\n")
