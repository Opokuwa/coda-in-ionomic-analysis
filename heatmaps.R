# =============================================================================
# Compositional Methods in Ionomic Data Analysis — Regression Significance Heatmaps
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Visualizes coefficient-level significance and beta direction from the
#   element-by-element linear models. Produces four heatmaps used in
#   Section 3.4 of the report:
#     - Raw proportion model significance
#     - CLR model significance
#     - Side-by-side significance comparison (Raw vs CLR)
#     - Beta direction + significance overlay (Raw vs CLR)
#
# PREREQUISITE: Run linreg.R first — requires raw_results, clr_results,
#               and elem_clean from that script.
#
# OUTPUT:
#   - Ionome_heatmaps_per_element.png      : coefficient significance (raw)
#   - Ionome_heatmaps_per_element_CLR.png  : coefficient significance (CLR)
#   - Ionome_heatmaps_raw_vs_clr.png       : side-by-side comparison
#   - Ionome_beta_direction_significance.png : beta direction with significance
# =============================================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(broom)
library(grid)


# =============================================================================
# 1. THEME
# =============================================================================

theme_thesis <- theme_bw(base_size = 16) +
  theme(
    plot.title  = element_text(size = 18, face = "bold"),
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13),
    legend.key.size = unit(0.7, "cm"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )


# =============================================================================
# 2. TIDY MODEL OUTPUTS
# =============================================================================
# Extract coefficient-level results from all element models into tidy format

tidy_from_list <- function(res_list, model_name) {
  bind_rows(lapply(res_list, function(x) {
    tidy(x$full) %>%
      mutate(
        term_clean = str_replace_all(term, "`", ""),
        model      = model_name,
        element    = x$element
      )
  }))
}

tidy_raw_all <- tidy_from_list(raw_results, "Raw")
tidy_clr_all <- tidy_from_list(clr_results, "CLR")


# =============================================================================
# 3. SIGNIFICANCE CATEGORIES
# =============================================================================
# Classify each predictor-element pair by p-value threshold.
# Region terms are excluded from the heatmaps (too many levels for readability).

categorize_sig <- function(tidy_df) {
  tidy_df %>%
    mutate(
      sig_cat = case_when(
        p.value < 0.001 ~ "p < 0.001",
        p.value < 0.01  ~ "p < 0.01",
        p.value < 0.05  ~ "p < 0.05",
        TRUE            ~ "ns"
      )
    ) %>%
    filter(term_clean != "(Intercept)",
           !str_detect(term_clean, "^Region")) %>%
    dplyr::select(response = element, term = term_clean, sig_cat)
}

tidy_raw <- categorize_sig(tidy_raw_all)
tidy_clr <- categorize_sig(tidy_clr_all)


# =============================================================================
# 4. HEATMAP: RAW PROPORTION SIGNIFICANCE
# =============================================================================

p1 <- ggplot(tidy_raw, aes(x = term, y = response, fill = sig_cat)) +
  geom_tile(color = "grey80") +
  scale_fill_manual(
    values = c("ns" = "#d9d9d9", "p < 0.05" = "#9ecae1",
               "p < 0.01" = "#4292c6", "p < 0.001" = "#084594"),
    name = "Significance"
  ) +
  labs(x = "Predictor", y = "Element (Raw proportions)",
       title = "Significant Effects: Raw Ionome Models") +
  theme_thesis


# =============================================================================
# 5. HEATMAP: CLR SIGNIFICANCE
# =============================================================================

p2 <- ggplot(tidy_clr, aes(x = term, y = response, fill = sig_cat)) +
  geom_tile(color = "grey80") +
  scale_fill_manual(
    values = c("ns" = "#d9d9d9", "p < 0.05" = "#fcae91",
               "p < 0.01" = "#fb6a4a", "p < 0.001" = "#cb181d"),
    name = "Significance"
  ) +
  labs(x = "Predictor", y = "Element (CLR coordinates)",
       title = "Significant Effects: CLR Ionome Models") +
  theme_thesis


# =============================================================================
# 6. HEATMAP: SIDE-BY-SIDE COMPARISON (RAW vs CLR)
# =============================================================================
# Left tile = Raw, Right tile = CLR for each predictor-element cell

tidy_raw$model <- "Raw"
tidy_clr$model <- "CLR"
coef_both <- bind_rows(tidy_raw, tidy_clr)
coef_both$model <- factor(coef_both$model, levels = c("Raw", "CLR"))

p3 <- ggplot(coef_both, aes(x = term, y = response)) +
  geom_tile(data = subset(coef_both, model == "Raw"),
            aes(fill = sig_cat), color = "black",
            width = 0.45, position = position_nudge(x = -0.25)) +
  geom_tile(data = subset(coef_both, model == "CLR"),
            aes(fill = sig_cat), color = "white",
            width = 0.45, position = position_nudge(x =  0.25)) +
  scale_fill_manual(
    values = c("ns" = "#d9d9d9", "p < 0.05" = "#fcae91",
               "p < 0.01" = "#fb6a4a", "p < 0.001" = "#cb181d"),
    name = "Significance"
  ) +
  labs(x = "Predictor (left = Raw, right = CLR)", y = "Element",
       title = "Environmental Effects on Ionome: Raw vs CLR Models") +
  theme_thesis


# =============================================================================
# 7. HEATMAP: BETA DIRECTION + SIGNIFICANCE
# =============================================================================
# Tile color shows direction of effect (positive/negative);
# asterisk (*) marks statistically significant cells (p < 0.05).

betas_raw <- tidy_raw_all %>%
  filter(term_clean != "(Intercept)") %>%
  transmute(response = element, term = term_clean, estimate, p.value, model = "Raw")

betas_clr <- tidy_clr_all %>%
  filter(term_clean != "(Intercept)") %>%
  transmute(response = element, term = term_clean, estimate, p.value, model = "CLR")

betas_both <- bind_rows(betas_raw, betas_clr) %>%
  mutate(
    sign_dir  = case_when(estimate > 0 ~ "positive",
                          estimate < 0 ~ "negative",
                          TRUE         ~ "zero"),
    sig       = p.value < 0.05,
    sig_label = ifelse(sig, "*", "")
  )
betas_both$model <- factor(betas_both$model, levels = c("Raw", "CLR"))

p4 <- ggplot(betas_both, aes(x = term, y = response)) +
  geom_tile(data = subset(betas_both, model == "Raw"),
            aes(fill = sign_dir), width = 0.45,
            position = position_nudge(x = -0.25), color = "grey80") +
  geom_tile(data = subset(betas_both, model == "CLR"),
            aes(fill = sign_dir), width = 0.45,
            position = position_nudge(x =  0.25), color = "grey80") +
  scale_fill_manual(
    values = c("negative" = "#e34a33", "zero" = "#f7f7f7", "positive" = "#3182bd"),
    name = "Beta direction"
  ) +
  geom_text(data = subset(betas_both, sig),
            aes(label = sig_label, x = term, y = response, group = model),
            position = position_nudge(
              x = ifelse(betas_both$model[betas_both$sig] == "Raw", -0.25, 0.25)
            ),
            size = 5, fontface = "bold", color = "white") +
  labs(x = "Predictor (left = Raw, right = CLR)", y = "Element",
       title = "Beta Direction & Significance: Raw vs CLR Ionome Models") +
  theme_thesis


# =============================================================================
# 8. EXPORT
# =============================================================================

ggsave("Ionome_heatmaps_per_element.png",          p1, width = 12, height = 8, dpi = 300)
ggsave("Ionome_heatmaps_per_element_CLR.png",      p2, width = 12, height = 8, dpi = 300)
ggsave("Ionome_heatmaps_raw_vs_clr.png",           p3, width = 12, height = 8, dpi = 300)
ggsave("Ionome_beta_direction_significance.png",   p4, width = 12, height = 8, dpi = 300)

cat("Saved 4 heatmap PNGs.\n")
