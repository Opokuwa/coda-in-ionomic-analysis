# =============================================================================
# Compositional Methods in Ionomic Data Analysis — Main Simulation Study
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Reproduces the main simulation study comparing Raw proportions,
#   log-transformed raw proportions (LogRaw), and CLR on Type I error
#   and statistical power. Self-contained — does not require linreg.R.
#
# REPRODUCES:
#   - Table 4.1  : Type I error and power under null and effect scenarios
#   - Figure 4.1 : Type I error rate by method and predictor type
#   - Figure 4.2 : Power by method and predictor type
#
# INPUT:   Amphipod dataset (used only to set realistic n and D values and
#          to estimate the baseline composition mu_hat)
# OUTPUT:
#   - table_4_1_simulation_results.csv
#   - figure_4_1_type1_error.png
#   - figure_4_2_power.png
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)

set.seed(123)


# =============================================================================
# 1. LOAD REAL DATA — EXTRACT SAMPLE SIZE AND BASELINE COMPOSITION
# =============================================================================
# The real amphipod dataset is used here only to:
#   (a) set n and D to match the empirical study
#   (b) estimate a realistic baseline mean composition (mu_hat) for simulation

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

X_raw <- as.matrix(dat[, elem_clean])
n <- nrow(X_raw)   # 113 — sample size for simulation
D <- ncol(X_raw)   # 28  — number of elements

# Convert raw concentrations to proportions; estimate baseline mean composition
X_prop  <- X_raw / rowSums(X_raw, na.rm = TRUE)
mu_hat  <- colMeans(X_prop, na.rm = TRUE)
mu_hat  <- pmax(mu_hat, 1e-8)   # avoid log(0) issues
mu_hat  <- mu_hat / sum(mu_hat)
names(mu_hat) <- elem_clean


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

# Dirichlet sampler via normalized gamma draws
rdirichlet_simple <- function(n, alpha) {
  out <- matrix(NA, nrow = n, ncol = length(alpha))
  for (i in seq_len(n)) {
    z       <- rgamma(length(alpha), shape = alpha, rate = 1)
    out[i,] <- z / sum(z)
  }
  out
}

# Simulate a compositional dataset under a log-linear model.
# The predictor shifts the log-scale mean of affected elements by beta_vec.
# Compositions are sampled from a Dirichlet distribution around the shifted mean.
simulate_composition <- function(n, mu, predictor_type = "binary",
                                  beta_vec = NULL, conc = 200,
                                  affected = 1:4) {
  D <- length(mu)
  x <- if (predictor_type == "binary") rbinom(n, 1, 0.5) else rnorm(n, 0, 1)
  if (is.null(beta_vec)) beta_vec <- rep(0, length(affected))

  Y <- matrix(NA, nrow = n, ncol = D)
  for (i in seq_len(n)) {
    effect        <- rep(0, D)
    effect[affected] <- beta_vec
    eta  <- log(mu) + x[i] * effect
    mu_i <- exp(eta); mu_i <- mu_i / sum(mu_i)
    Y[i,] <- rdirichlet_simple(1, conc * mu_i)
  }
  colnames(Y) <- names(mu)
  list(Y = Y, x = x)
}

# Analysis functions: return element-wise p-values for each method
get_raw_pvals <- function(Y, x) {
  apply(Y, 2, function(y) anova(lm(y ~ 1), lm(y ~ x))$`Pr(>F)`[2])
}

get_lograw_pvals <- function(Y, x) {
  Y_log <- log(pmax(Y, 1e-8))
  apply(Y_log, 2, function(y) anova(lm(y ~ 1), lm(y ~ x))$`Pr(>F)`[2])
}

get_clr_pvals <- function(Y, x) {
  Y_safe <- pmax(Y, 1e-8)
  Y_clr  <- log(Y_safe) - rowMeans(log(Y_safe))
  apply(Y_clr, 2, function(y) anova(lm(y ~ 1), lm(y ~ x))$`Pr(>F)`[2])
}


# =============================================================================
# 3. SIMULATION DRIVER
# =============================================================================
# Runs nsim iterations under a given scenario, computing element-wise rejection
# rates for Raw, LogRaw, and CLR — both unadjusted and Bonferroni-adjusted.
# Returns a summary data frame suitable for Table 4.1.

run_sim_study <- function(nsim = 1000, n, predictor_type = "binary",
                           beta_scale = 0, conc = 200, alpha_level = 0.05,
                           affected = 1:4, D) {
  methods <- c("Raw", "LogRaw", "CLR")
  raw_p <- lograw_p <- clr_p <-
    raw_p_adj <- lograw_p_adj <- clr_p_adj <- matrix(NA, nrow = nsim, ncol = D)
  colnames(raw_p) <- colnames(lograw_p) <- colnames(clr_p) <-
    colnames(raw_p_adj) <- colnames(lograw_p_adj) <- colnames(clr_p_adj) <- elem_clean

  for (s in seq_len(nsim)) {
    beta_vec <- if (beta_scale == 0) {
      rep(0, length(affected))
    } else {
      sample(c(-1, 1), length(affected), replace = TRUE) *
        runif(length(affected), 0.2, 1.0) * beta_scale
    }

    sim <- simulate_composition(n, mu_hat, predictor_type, beta_vec, conc, affected)

    raw_p[s,]    <- get_raw_pvals(sim$Y, sim$x)
    lograw_p[s,] <- get_lograw_pvals(sim$Y, sim$x)
    clr_p[s,]    <- get_clr_pvals(sim$Y, sim$x)

    raw_p_adj[s,]    <- p.adjust(raw_p[s,],    "bonferroni")
    lograw_p_adj[s,] <- p.adjust(lograw_p[s,], "bonferroni")
    clr_p_adj[s,]    <- p.adjust(clr_p[s,],    "bonferroni")
  }

  raw_rate      <- colMeans(raw_p    < alpha_level, na.rm = TRUE)
  lograw_rate   <- colMeans(lograw_p < alpha_level, na.rm = TRUE)
  clr_rate      <- colMeans(clr_p    < alpha_level, na.rm = TRUE)
  raw_bonf      <- colMeans(raw_p_adj    < alpha_level, na.rm = TRUE)
  lograw_bonf   <- colMeans(lograw_p_adj < alpha_level, na.rm = TRUE)
  clr_bonf      <- colMeans(clr_p_adj    < alpha_level, na.rm = TRUE)

  unaffected <- setdiff(seq_len(D), affected)
  is_null    <- beta_scale == 0

  data.frame(
    predictor  = predictor_type,
    scenario   = ifelse(is_null, "null", "effect"),
    beta_scale = beta_scale,
    adjustment = rep(c("unadjusted", "bonferroni"), each = 3),
    method     = rep(methods, 2),
    power      = if (is_null) rep(NA_real_, 6) else
      c(mean(raw_rate[affected]),   mean(lograw_rate[affected]),   mean(clr_rate[affected]),
        mean(raw_bonf[affected]),   mean(lograw_bonf[affected]),   mean(clr_bonf[affected])),
    type1      = c(mean(raw_rate[if (is_null) seq_len(D) else unaffected]),
                   mean(lograw_rate[if (is_null) seq_len(D) else unaffected]),
                   mean(clr_rate[if (is_null) seq_len(D) else unaffected]),
                   mean(raw_bonf[if (is_null) seq_len(D) else unaffected]),
                   mean(lograw_bonf[if (is_null) seq_len(D) else unaffected]),
                   mean(clr_bonf[if (is_null) seq_len(D) else unaffected]))
  )
}


# =============================================================================
# 4. RUN ALL SCENARIOS
# =============================================================================
# Four scenarios: null/effect × binary/continuous predictor
# 1,000 simulations each; concentration parameter fixed at 200

set.seed(123)

null_binary     <- run_sim_study(1000, n, "binary",     0,    200, 0.05, 1:4, D)
null_continuous <- run_sim_study(1000, n, "continuous", 0,    200, 0.05, 1:4, D)
eff_binary      <- run_sim_study(1000, n, "binary",    -0.8,  200, 0.05, 1:4, D)
eff_continuous  <- run_sim_study(1000, n, "continuous",-0.8,  200, 0.05, 1:4, D)

sim_summary <- bind_rows(null_binary, null_continuous, eff_binary, eff_continuous)


# =============================================================================
# 5. TABLE 4.1
# =============================================================================

table_4_1 <- sim_summary %>%
  mutate(predictor  = ifelse(predictor == "binary", "Binary", "Continuous"),
         adjustment = ifelse(adjustment == "unadjusted", "Unadjusted", "Bonferroni"))

print(table_4_1)
write.csv(table_4_1, "table_4_1_simulation_results.csv", row.names = FALSE)
cat("Table 4.1 saved to table_4_1_simulation_results.csv\n")


# =============================================================================
# 6. FIGURE 4.1 — TYPE I ERROR
# =============================================================================

sim_summary$method     <- factor(sim_summary$method,     levels = c("Raw", "LogRaw", "CLR"))
sim_summary$predictor  <- factor(sim_summary$predictor,  levels = c("binary", "continuous"))
sim_summary$adjustment <- factor(sim_summary$adjustment, levels = c("unadjusted", "bonferroni"))

fig_4_1 <- ggplot(sim_summary %>% filter(scenario == "null"),
                  aes(x = method, y = type1, fill = method)) +
  geom_col(width = 0.7) +
  facet_grid(adjustment ~ predictor) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(title = "Type I Error by Method and Predictor Type",
       x = "Method", y = "Type I error rate") +
  theme_bw() + theme(legend.position = "none")

print(fig_4_1)
ggsave("figure_4_1_type1_error.png", fig_4_1, width = 10, height = 7, dpi = 300)
cat("Figure 4.1 saved to figure_4_1_type1_error.png\n")


# =============================================================================
# 7. FIGURE 4.2 — POWER
# =============================================================================

fig_4_2 <- ggplot(sim_summary %>% filter(scenario == "effect"),
                  aes(x = method, y = power, fill = method)) +
  geom_col(width = 0.7) +
  facet_grid(adjustment ~ predictor) +
  labs(title = "Power by Method and Predictor Type",
       x = "Method", y = "Power") +
  theme_bw() + theme(legend.position = "none")

print(fig_4_2)
ggsave("figure_4_2_power.png", fig_4_2, width = 10, height = 7, dpi = 300)
cat("Figure 4.2 saved to figure_4_2_power.png\n")
