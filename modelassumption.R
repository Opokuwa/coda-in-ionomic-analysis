# =============================================================================
# Compositional Methods in Ionomic Data Analysis — Model Assumption Checks
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Validates linear model assumptions for raw and CLR models across all 28
#   elements. Produces diagnostic plots (residuals, Q-Q, Cook's distance) and
#   a summary table of Shapiro-Wilk and Durbin-Watson test results.
#
# PREREQUISITE: Run linreg.R first — requires raw_results, clr_results,
#               and elem_clean from that script.
#
# OUTPUT:
#   - assumption_diagnostics.pdf  : full diagnostic plots for all elements
#   - diagnostics_RAW_<el>.png    : individual PNG per element (raw model)
#   - diagnostics_CLR_<el>.png    : individual PNG per element (CLR model)
#   - assumption_summary_df       : data frame of Shapiro-Wilk + DW p-values
# =============================================================================

library(lmtest)


# =============================================================================
# 1. DIAGNOSTIC PLOT FUNCTION
# =============================================================================
# Produces a 2×2 panel: Residuals vs Fitted, Normal Q-Q, Residuals vs Order,
# Cook's Distance. The Cook's distance threshold line is drawn at 4/n.

check_lm_assumptions <- function(fit, main_title = "") {
  op <- par(mfrow = c(2, 2))

  # Residuals vs Fitted — checks linearity and homoscedasticity
  plot(fitted(fit), resid(fit),
       xlab = "Fitted values",
       ylab = "Residuals",
       main = paste("Residuals vs Fitted", main_title))
  abline(h = 0, col = "red")

  # Normal Q-Q — checks normality of residuals
  qqnorm(resid(fit), main = paste("Normal Q-Q", main_title))
  qqline(resid(fit), col = "red")

  # Residuals vs Observation Order — checks independence
  plot(resid(fit), type = "b",
       xlab = "Observation order",
       ylab = "Residuals",
       main = paste("Residuals vs Order", main_title))
  abline(h = 0, col = "red")

  # Cook's Distance — identifies influential observations
  plot(cooks.distance(fit),
       ylab = "Cook's distance",
       xlab = "Observation",
       main = paste("Cook's Distance", main_title),
       type = "h")
  abline(h = 4 / length(resid(fit)), col = "red", lty = 2)

  par(op)
}


# =============================================================================
# 2. SINGLE-ELEMENT EXAMPLE
# =============================================================================
# Quick visual check before running the full batch

fit_raw_ca <- raw_results[[which(elem_clean == "Ca")]]$full
check_lm_assumptions(fit_raw_ca, main_title = "RAW Ca")

fit_clr_ca <- clr_results[[which(elem_clean == "Ca")]]$full
check_lm_assumptions(fit_clr_ca, main_title = "CLR Ca")


# =============================================================================
# 3. BATCH PDF OUTPUT — ALL ELEMENTS
# =============================================================================
# Saves all 28×2 diagnostic panels (raw + CLR) to a single PDF

pdf("assumption_diagnostics.pdf")
for (el in elem_clean) {
  fit_raw <- raw_results[[which(elem_clean == el)]]$full
  check_lm_assumptions(fit_raw, main_title = paste("RAW", el))

  fit_clr <- clr_results[[which(elem_clean == el)]]$full
  check_lm_assumptions(fit_clr, main_title = paste("CLR", el))
}
dev.off()
cat("Saved diagnostic plots to assumption_diagnostics.pdf\n")


# =============================================================================
# 4. INDIVIDUAL PNG OUTPUTS — ALL ELEMENTS
# =============================================================================
# Exports one PNG per element per model for selective use in the report

par(mfrow = c(2, 2), cex = 1.2)

for (el in elem_clean) {

  png(filename = paste0("diagnostics_RAW_", el, ".png"),
      width = 1200, height = 1000, res = 150)
  fit_raw <- raw_results[[which(elem_clean == el)]]$full
  check_lm_assumptions(fit_raw, main_title = paste("RAW", el))
  dev.off()

  png(filename = paste0("diagnostics_CLR_", el, ".png"),
      width = 1200, height = 1000, res = 150)
  fit_clr <- clr_results[[which(elem_clean == el)]]$full
  check_lm_assumptions(fit_clr, main_title = paste("CLR", el))
  dev.off()
}

cat("Saved individual diagnostic PNGs for all elements.\n")


# =============================================================================
# 5. FORMAL TEST SUMMARY TABLE
# =============================================================================
# Shapiro-Wilk (normality) and Durbin-Watson (independence) p-values
# for both raw and CLR models across all 28 elements.
# Sorted by smallest raw DW p-value to highlight problematic cases.

get_shapiro_p <- function(fit) shapiro.test(resid(fit))$p.value
get_dw_p      <- function(fit) lmtest::dwtest(fit)$p.value

assumption_summary <- lapply(seq_along(elem_clean), function(i) {
  el      <- elem_clean[i]
  raw_fit <- raw_results[[i]]$full
  clr_fit <- clr_results[[i]]$full

  data.frame(
    element       = el,
    raw_shapiro_p = get_shapiro_p(raw_fit),
    raw_dw_p      = get_dw_p(raw_fit),
    clr_shapiro_p = get_shapiro_p(clr_fit),
    clr_dw_p      = get_dw_p(clr_fit)
  )
})

assumption_summary_df <- do.call(rbind, assumption_summary)
assumption_summary_df <- assumption_summary_df[order(assumption_summary_df$raw_dw_p), ]

print(assumption_summary_df)
