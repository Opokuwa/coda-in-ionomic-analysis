# =============================================================================
# Compositional Methods in Ionomic Data Analysis — Compositional Bar Plots
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Visualizes the elemental composition of amphipod samples under two
#   representations — raw proportions and CLR-transformed coordinates —
#   grouped by Aquifer and Sex. Produces Figures 3.1–3.4 in the report.
#
# INPUT:   Amphipod ionomic dataset (113 samples × 28 elements)
# OUTPUT:
#   - barplot_raw_aquifer.png   : stacked bar — raw proportions by Aquifer
#   - barplot_raw_sex.png       : stacked bar — raw proportions by Sex
#   - barplot_clr_aquifer.png   : grouped bar — mean CLR by Aquifer
#   - barplot_clr_sex.png       : grouped bar — mean CLR by Sex
# =============================================================================

library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(scales)
library(compositions)
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

# Consistent color palette — one color per element
my_cols <- c(
  "Al" = "#e41a1c", "As" = "#32CD32", "B"  = "#0000FF", "Ba" = "#984ea3",
  "Be" = "#ff7f00", "Bi" = "#ffff33", "Ca" = "#8B0000", "Cd" = "#fccde5",
  "Co" = "#1b9e77", "Cr" = "#E5E5E5", "Cu" = "#7570b3", "Fe" = "#EEAD0E",
  "K"  = "#66a61e", "Li" = "#B3EE3A", "Mg" = "#a6761d", "Mn" = "#666666",
  "Mo" = "#8dd3c7", "Na" = "#FF1493", "Ni" = "#80b1d3", "P"  = "#7A67EE",
  "Pb" = "#b3de69", "S"  = "#2F4F4F", "Se" = "#00CED1", "Si" = "#bc80bd",
  "Sr" = "#ccebc5", "Tl" = "#ffed6f", "V"  = "#1f78b4", "Zn" = "#FF00FF"
)

# Global ggplot theme
theme_set(
  theme_bw(base_size = 18) +
    theme(
      plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title    = element_text(size = 18, face = "bold"),
      axis.text     = element_text(size = 15, color = "black"),
      legend.title  = element_text(size = 16, face = "bold"),
      legend.text   = element_text(size = 13, color = "black"),
      legend.key.size = unit(0.9, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text    = element_text(size = 16, face = "bold")
    )
)


# =============================================================================
# 2. RAW PROPORTION BAR PLOTS
# =============================================================================
# Stacked bar charts showing mean elemental proportions by Aquifer and Sex.
# Each bar sums to 1; colors distinguish elements.

dat_rel <- dat %>%
  filter(!is.na(Sex)) %>%
  rowwise() %>%
  mutate(Total_ions = sum(c_across(all_of(elem_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_longer(cols = all_of(elem_cols),
               names_to  = "Element",
               values_to = "Conc") %>%
  filter(!is.na(Conc), Total_ions > 0) %>%
  mutate(
    Prop    = Conc / Total_ions,
    Element = gsub(" \\(ug/g\\)", "", Element)
  )

# Mean proportion per Aquifer × Element
aquifer_comp <- dat_rel %>%
  group_by(Aquifer, Element) %>%
  summarise(Mean_prop = mean(Prop, na.rm = TRUE), .groups = "drop")

# Mean proportion per Sex × Element
gender_comp <- dat_rel %>%
  group_by(Sex, Element) %>%
  summarise(Mean_prop = mean(Prop, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(aquifer_comp, aes(x = Aquifer, y = Mean_prop, fill = Element)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "Average Compositions by Aquifer",
       x = "Aquifer", y = "Proportion", fill = "Element") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 2))

p2 <- ggplot(gender_comp, aes(x = Sex, y = Mean_prop, fill = Element)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = my_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "Average Compositions by Sex",
       x = "Sex", y = "Proportion", fill = "Element") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 2))

p1
p2


# =============================================================================
# 3. CLR BAR PLOTS
# =============================================================================
# Grouped bar charts of mean CLR-transformed values by Aquifer and Sex.
# Unlike raw proportions, CLR values are not bounded to [0,1] and can be
# negative; the dashed line at 0 marks the geometric mean of the composition.

dat_sub <- dat %>% filter(!is.na(Aquifer), !is.na(Sex))
X_raw   <- dat_sub %>% dplyr::select(all_of(elem_cols))
X_comp  <- acomp(X_raw)
X_clr   <- clr(X_comp)

clr_df <- as.data.frame(X_clr)
names(clr_df) <- gsub(" \\(ug/g\\)", "", elem_cols)
clr_df <- clr_df %>%
  mutate(Aquifer = dat_sub$Aquifer,
         Sex     = dat_sub$Sex)

elem_clean_local <- gsub(" \\(ug/g\\)", "", elem_cols)

# CLR by Aquifer
aquifer_clr <- clr_df %>%
  pivot_longer(cols = all_of(elem_clean_local),
               names_to = "Element", values_to = "CLR") %>%
  group_by(Aquifer, Element) %>%
  summarise(Mean_CLR = mean(CLR, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(aquifer_clr, aes(x = Element, y = Mean_CLR, fill = Aquifer)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Mean CLR-Transformed Ionome by Aquifer",
       x = "Element", y = "Mean CLR", fill = "Aquifer") +
  coord_flip() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(0.6, "cm")) +
  guides(fill = guide_legend(nrow = 1))

# CLR by Sex
sex_clr <- clr_df %>%
  pivot_longer(cols = all_of(elem_clean_local),
               names_to = "Element", values_to = "CLR") %>%
  group_by(Sex, Element) %>%
  summarise(Mean_CLR = mean(CLR, na.rm = TRUE), .groups = "drop")

p4 <- ggplot(sex_clr, aes(x = Element, y = Mean_CLR, fill = Sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Mean CLR-Transformed Ionome by Sex",
       x = "Element", y = "Mean CLR", fill = "Sex") +
  coord_flip() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(0.6, "cm")) +
  guides(fill = guide_legend(nrow = 1))

p3
p4


# =============================================================================
# 4. EXPORT
# =============================================================================

ggsave("barplot_raw_aquifer.png", p1, width = 10, height = 8, dpi = 300)
ggsave("barplot_raw_sex.png",     p2, width = 10, height = 8, dpi = 300)
ggsave("barplot_clr_aquifer.png", p3, width = 10, height = 8, dpi = 300)
ggsave("barplot_clr_sex.png",     p4, width = 10, height = 8, dpi = 300)

cat("Saved 4 bar plot PNGs.\n")
