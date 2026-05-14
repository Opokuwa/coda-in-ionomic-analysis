# =============================================================================
# Compositional Methods in Ionomic Data Analysis — ANOVA-Level Heatmaps
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Visualizes overall predictor significance (ANOVA F-tests) rather than
#   individual coefficient significance. Shows whether each predictor
#   (Sex, Aquifer, Region) has a significant *overall* effect on each element
#   under raw and CLR models.
#
# PREREQUISITE: Run linreg.R first — requires raw_results, clr_results,
#               and elem_clean from that script.
#
# OUTPUT:
#   - Ionome_heatmaps_overalleffect.pdf — three ANOVA heatmaps
# =============================================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(grid)


# =============================================================================
# 1. THEME
# =============================================================================

theme_set(
  theme_bw(base_size = 18) +
    theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title    = element_text(size = 18, face = "bold"),
      axis.text     = element_text(size = 14, color = "black"),
      legend.title  = element_text(size = 14, face = "bold"),
      legend.text   = element_text(size = 12, color = "black"),
      legend.key.size = unit(0.8, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
)


# =============================================================================
# 2. EXTRACT ANOVA P-VALUES
# =============================================================================
# Each row is one element; columns are the overall ANOVA p-values for
# Sex, Aquifer, and Region (from the likelihood ratio tests in linreg.R).

get_anova_p <- function(anova_obj) as.data.frame(anova_obj)$`Pr(>F)`[2]

anova_raw_df <- map_dfr(raw_results, function(x) {
  tibble(element = x$element,
         Sex     = get_anova_p(x$anova_no_sex),
         Aquifer = get_anova_p(x$anova_no_aq),
         Region  = get_anova_p(x$anova_no_reg))
})

anova_clr_df <- map_dfr(clr_results, function(x) {
  tibble(element = x$element,
         Sex     = get_anova_p(x$anova_no_sex),
         Aquifer = get_anova_p(x$anova_no_aq),
         Region  = get_anova_p(x$anova_no_reg))
})


# =============================================================================
# 3. PIVOT TO LONG FORMAT AND CATEGORIZE
# =============================================================================

categorize_anova <- function(df) {
  df %>%
    pivot_longer(cols = c(Sex, Aquifer, Region),
                 names_to  = "predictor",
                 values_to = "p_value") %>%
    mutate(
      sig_cat = case_when(
        p_value < 0.001 ~ "p < 0.001",
        p_value < 0.01  ~ "p < 0.01",
        p_value < 0.05  ~ "p < 0.05",
        TRUE            ~ "ns"
      ),
      sig_cat = factor(sig_cat, levels = c("ns", "p < 0.05", "p < 0.01", "p < 0.001"))
    )
}

anova_raw_long <- categorize_anova(anova_raw_df)
anova_clr_long <- categorize_anova(anova_clr_df)

# Preserve original element order along the y-axis
anova_raw_long$element <- factor(anova_raw_long$element,
                                  levels = rev(unique(anova_raw_long$element)))
anova_clr_long$element <- factor(anova_clr_long$element,
                                  levels = rev(unique(anova_clr_long$element)))


# =============================================================================
# 4. HEATMAP: RAW PROPORTION MODELS
# =============================================================================

p5 <- ggplot(anova_raw_long, aes(x = predictor, y = element, fill = sig_cat)) +
  geom_tile(color = "grey70", linewidth = 0.6) +
  scale_fill_manual(
    values = c("ns" = "#d9d9d9", "p < 0.05" = "#9ecae1",
               "p < 0.01" = "#4292c6", "p < 0.001" = "#084594"),
    name = "Overall p-value",
    drop = FALSE
  ) +
  labs(x = "Predictor", y = "Element",
       title = "Overall ANOVA p-values: Raw Ionome Models") +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 13),
        legend.position = "right")


# =============================================================================
# 5. HEATMAP: CLR MODELS
# =============================================================================

p6 <- ggplot(anova_clr_long, aes(x = predictor, y = element, fill = sig_cat)) +
  geom_tile(color = "grey70", linewidth = 0.6) +
  scale_fill_manual(
    values = c("ns" = "#d9d9d9", "p < 0.05" = "#fcae91",
               "p < 0.01" = "#fb6a4a", "p < 0.001" = "#cb181d"),
    name = "Overall p-value",
    drop = FALSE
  ) +
  labs(x = "Predictor", y = "Element",
       title = "Overall ANOVA p-values: CLR Ionome Models") +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 13),
        legend.position = "right")


# =============================================================================
# 6. HEATMAP: RAW vs CLR COMPARISON
# =============================================================================
# Split tiles — left = Raw, right = CLR — for direct cell-by-cell comparison

anova_raw_long$model <- "Raw"
anova_clr_long$model <- "CLR"
anova_both_long <- bind_rows(anova_raw_long, anova_clr_long)
anova_both_long$model   <- factor(anova_both_long$model,   levels = c("Raw", "CLR"))
anova_both_long$element <- factor(anova_both_long$element,
                                   levels = rev(unique(anova_both_long$element)))

p7 <- ggplot(anova_both_long, aes(x = predictor, y = element)) +
  geom_tile(data = subset(anova_both_long, model == "Raw"),
            aes(fill = sig_cat), width = 0.42, height = 0.85,
            position = position_nudge(x = -0.22),
            color = "black", linewidth = 0.5) +
  geom_tile(data = subset(anova_both_long, model == "CLR"),
            aes(fill = sig_cat), width = 0.42, height = 0.85,
            position = position_nudge(x =  0.22),
            color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("ns" = "#d9d9d9", "p < 0.05" = "#fcae91",
               "p < 0.01" = "#fb6a4a", "p < 0.001" = "#cb181d"),
    name = "Overall p-value",
    drop = FALSE
  ) +
  labs(x = "Predictor", y = "Element",
       title   = "Overall ANOVA p-values by Predictor: Raw vs CLR Models",
       subtitle = "Left tile = Raw  |  Right tile = CLR") +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 13),
        legend.position = "right",
        plot.subtitle = element_text(size = 14, hjust = 0.5))


# =============================================================================
# 7. EXPORT
# =============================================================================

pdf("Ionome_heatmaps_overalleffect.pdf", width = 12, height = 8)
print(p5)
print(p6)
print(p7)
dev.off()

cat("Saved 3 ANOVA heatmaps to Ionome_heatmaps_overalleffect.pdf\n")
