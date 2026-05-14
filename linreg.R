# =============================================================================
# Compositional Methods in Ionomic Data Analysis — Linear Regression Models
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Fits element-by-element linear models under two representations:
#   (1) Raw proportions and (2) CLR-transformed coordinates.
#   Compares predictor significance (Sex, Aquifer, Region) across both
#   approaches as the primary real-data analysis in the report.
#
# INPUT:   Amphipod ionomic dataset (113 samples × 28 elements)
# OUTPUT:  raw_results, clr_results — lists of fitted models used downstream
#
# RUN ORDER: Run this script first. modelassumption.R, heatmaps.R, and
#            anovaheatmaps.R all depend on objects created here.
# =============================================================================

library(compositions)
library(dplyr)
library(broom)
library(readxl)

# =============================================================================
# 1. DATA LOADING AND PREPARATION
# =============================================================================

# Update the path below to point to your local copy of the dataset
fn  <- "OK Springs master datasheet.xlsx"
dat <- read_excel(fn, sheet = "Pods ICP")

# Element column names as they appear in the raw data
elem_cols <- c("Al (ug/g)", "As (ug/g)", "B (ug/g)",  "Ba (ug/g)", "Be (ug/g)",
               "Bi (ug/g)", "Ca (ug/g)", "Cd (ug/g)", "Co (ug/g)", "Cr (ug/g)",
               "Cu (ug/g)", "Fe (ug/g)", "K (ug/g)",  "Li (ug/g)", "Mg (ug/g)",
               "Mn (ug/g)", "Mo (ug/g)", "Na (ug/g)", "Ni (ug/g)", "P (ug/g)",
               "Pb (ug/g)", "S (ug/g)",  "Se (ug/g)", "Si (ug/g)", "Sr (ug/g)",
               "Tl (ug/g)", "V (ug/g)",  "Zn (ug/g)")

# Rename columns to clean element symbols (e.g., "Al (ug/g)" -> "Al")
elem_clean <- gsub(" \\(ug/g\\)", "", elem_cols)
names(dat)[match(elem_cols, names(dat))] <- elem_clean


# =============================================================================
# 2. REPRESENTATION MATRICES
# =============================================================================

# Raw elemental concentrations
X_raw  <- dat[, elem_clean]

# Raw proportions (row-wise closure to sum = 1)
Y_prop <- as.matrix(X_raw / rowSums(X_raw))

# CLR-transformed proportions (Aitchison geometry)
# CLR(x) = log(x) - mean(log(x)) — maps compositions to real space
X_clr  <- clr(acomp(Y_prop))
Y_clr  <- as.matrix(X_clr)

# Biological predictors
X_pred <- dat %>%
  dplyr::select(Sex, Aquifer, Region) %>%
  mutate(
    Sex     = factor(Sex),
    Aquifer = factor(Aquifer),
    Region  = factor(Region)
  )


# =============================================================================
# 3. ELEMENT-BY-ELEMENT LINEAR MODELS: RAW PROPORTIONS
# =============================================================================
# For each element, fit a full model (Sex + Aquifer + Region) and three
# reduced models to enable likelihood ratio tests for each predictor.

raw_results <- lapply(elem_clean, function(el) {
  full   <- lm(Y_prop[, el] ~ Sex + Aquifer + Region, data = X_pred)
  no_sex <- lm(Y_prop[, el] ~ Aquifer + Region,       data = X_pred)
  no_aq  <- lm(Y_prop[, el] ~ Sex + Region,           data = X_pred)
  no_reg <- lm(Y_prop[, el] ~ Sex + Aquifer,          data = X_pred)

  list(
    element      = el,
    full         = full,
    no_sex       = no_sex,
    no_aq        = no_aq,
    no_reg       = no_reg,
    anova_no_sex = anova(no_sex, full),
    anova_no_aq  = anova(no_aq,  full),
    anova_no_reg = anova(no_reg, full)
  )
})


# --- Example: Inspect results for specific elements -----------------------

# Zinc
raw_results[[which(elem_clean == "Zn")]]$anova_no_sex
raw_results[[which(elem_clean == "Zn")]]$anova_no_aq
raw_results[[which(elem_clean == "Zn")]]$anova_no_reg
summary(raw_results[[which(elem_clean == "Zn")]]$full)

# Calcium
raw_results[[which(elem_clean == "Ca")]]$anova_no_sex
raw_results[[which(elem_clean == "Ca")]]$anova_no_aq
raw_results[[which(elem_clean == "Ca")]]$anova_no_reg
summary(raw_results[[which(elem_clean == "Ca")]]$full)


# =============================================================================
# 4. ELEMENT-BY-ELEMENT LINEAR MODELS: CLR COORDINATES
# =============================================================================
# Identical structure to the raw models, applied to CLR-transformed values.
# CLR preserves the relative structure of the composition and maps it to
# an unconstrained real-valued space suitable for standard regression.

clr_results <- lapply(elem_clean, function(el) {
  full   <- lm(Y_clr[, el] ~ Sex + Aquifer + Region, data = X_pred)
  no_sex <- lm(Y_clr[, el] ~ Aquifer + Region,       data = X_pred)
  no_aq  <- lm(Y_clr[, el] ~ Sex + Region,           data = X_pred)
  no_reg <- lm(Y_clr[, el] ~ Sex + Aquifer,          data = X_pred)

  list(
    element      = el,
    full         = full,
    no_sex       = no_sex,
    no_aq        = no_aq,
    no_reg       = no_reg,
    anova_no_sex = anova(no_sex, full),
    anova_no_aq  = anova(no_aq,  full),
    anova_no_reg = anova(no_reg, full)
  )
})


# --- Example: Inspect CLR results for specific elements -------------------

# Zinc
clr_results[[which(elem_clean == "Zn")]]$anova_no_sex
clr_results[[which(elem_clean == "Zn")]]$anova_no_aq
clr_results[[which(elem_clean == "Zn")]]$anova_no_reg
summary(clr_results[[which(elem_clean == "Zn")]]$full)

# Calcium
clr_results[[which(elem_clean == "Ca")]]$anova_no_sex
clr_results[[which(elem_clean == "Ca")]]$anova_no_aq
clr_results[[which(elem_clean == "Ca")]]$anova_no_reg
summary(clr_results[[which(elem_clean == "Ca")]]$full)
