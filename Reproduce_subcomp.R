# =============================================================================
# Compositional Methods in Ionomic Data Analysis — Subcompositional Coherence
# Master's Report | Oklahoma State University, April 2026
# Author: Miriam Tenkorang | Supervisor: Dr. Pratyaydipta Rudra
# =============================================================================
# PURPOSE:
#   Demonstrates subcompositional coherence — the property that results should
#   remain consistent when a subset of components is analyzed rather than the
#   full composition. Shows that Raw and CLR violate this property while
#   pairwise log-ratios (PLR) preserve it.
#
#   Fully self-contained: no external data required.
#
# REPRODUCES:
#   - Table 4.2  : Estimated effect of binary predictor under Full, SubA, SubB
#                  compositions, for Raw, CLR, and PLR methods
#   - Figure 4.3 : Graphical display of subcompositional coherence
#
# OUTPUT:
#   - table_4_2_subcomp_results.csv
#   - figure_4_3_subcomp_coherence.png
# =============================================================================

library(ggplot2)
library(dplyr)

set.seed(42)


# =============================================================================
# 1. SIMULATION SETUP
# =============================================================================

n          <- 500          # number of observations
D          <- 25           # number of compositional components (C1, ..., C25)
kappa      <- 80           # Dirichlet concentration parameter
comp_names <- paste0("C", 1:D)

# Baseline log-scale parameters (intercept-only, equal for all components)
alpha <- rep(0, D)

# Effect vector: strong effect on C2–C20, weak effect on C1, none on C21–C25
# This structure creates two natural subcompositions for comparison:
#   SubA = {C1, C2, C3, C4, C5, C25}  — mixes strong and weak effects
#   SubB = {C1, C21, C22, C23, C24, C25}  — mixes weak and null effects
beta <- c(0.5, rep(3.0, 19), rep(0.0, 5))

# Binary predictor
X <- rbinom(n, size = 1, prob = 0.5)


# =============================================================================
# 2. GENERATE COMPOSITIONAL DATA
# =============================================================================

# Compute row-specific mean compositions under the log-linear model
eta_mat     <- outer(X, beta) + matrix(alpha, nrow = n, ncol = D, byrow = TRUE)
exp_eta_mat <- exp(eta_mat)
P           <- exp_eta_mat / rowSums(exp_eta_mat)

# Dirichlet sampler via normalized gamma draws
rdirichlet_simple <- function(alpha_vec) {
  z <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
  z / sum(z)
}

comp <- t(sapply(seq_len(n), function(i) rdirichlet_simple(kappa * P[i,])))
colnames(comp) <- comp_names

stopifnot(max(abs(rowSums(comp) - 1)) < 1e-9)


# =============================================================================
# 3. DEFINE COMPOSITIONAL CONTEXTS
# =============================================================================
# Each context represents a different subset of components that could be
# analyzed if a researcher focuses on a subgroup of elements.

idx_full <- 1:25                          # full 25-component composition
idx_subA <- c(1, 2, 3, 4, 5, 25)         # strong-effect-dominated subset
idx_subB <- c(1, 21, 22, 23, 24, 25)     # weak/null-effect subset


# =============================================================================
# 4. ANALYSIS TARGET FUNCTIONS
# =============================================================================
# For each method, compute the response value for component C1 after reclosing
# to the selected subset of components.

# Raw proportion of C1 within the selected composition
raw_prop_C1 <- function(mat, idx) {
  sub           <- mat[, idx, drop = FALSE]
  sub_reclosed  <- sub / rowSums(sub)
  sub_reclosed[, which(idx == 1)]
}

# CLR-transformed value of C1 within the selected composition
clr_C1 <- function(mat, idx) {
  sub          <- mat[, idx, drop = FALSE]
  sub_reclosed <- sub / rowSums(sub)
  log_sub      <- log(sub_reclosed)
  clr_sub      <- log_sub - rowMeans(log_sub)
  clr_sub[, which(idx == 1)]
}

# Pairwise log-ratio log(C1 / C25)
# This ratio is invariant to the choice of subcomposition — that is the key
# property being demonstrated. It remains the same regardless of which other
# components are included.
plr_C1_C25 <- function(mat, idx = NULL) {
  if (is.null(idx)) {
    return(log(mat[, 1]) - log(mat[, 25]))
  }
  sub         <- mat[, idx, drop = FALSE]
  sub_reclosed <- sub / rowSums(sub)
  pos1  <- which(idx == 1)
  pos25 <- which(idx == 25)
  log(sub_reclosed[, pos1]) - log(sub_reclosed[, pos25])
}


# =============================================================================
# 5. COMPUTE RESPONSE VARIABLES IN EACH CONTEXT
# =============================================================================

y_raw_full <- raw_prop_C1(comp, idx_full)
y_raw_subA <- raw_prop_C1(comp, idx_subA)
y_raw_subB <- raw_prop_C1(comp, idx_subB)

y_clr_full <- clr_C1(comp, idx_full)
y_clr_subA <- clr_C1(comp, idx_subA)
y_clr_subB <- clr_C1(comp, idx_subB)

y_plr_full <- plr_C1_C25(comp)
y_plr_subA <- plr_C1_C25(comp, idx_subA)
y_plr_subB <- plr_C1_C25(comp, idx_subB)

# Verify PLR invariance numerically (differences should be essentially zero)
cat("Max |PLR_full - PLR_subA|:", max(abs(y_plr_full - y_plr_subA)), "\n")
cat("Max |PLR_full - PLR_subB|:", max(abs(y_plr_full - y_plr_subB)), "\n\n")


# =============================================================================
# 6. FIT LINEAR MODELS
# =============================================================================
# The coefficient of X from each model is the estimated effect reported in Table 4.2.

m_raw_full <- lm(y_raw_full ~ X); m_raw_subA <- lm(y_raw_subA ~ X); m_raw_subB <- lm(y_raw_subB ~ X)
m_clr_full <- lm(y_clr_full ~ X); m_clr_subA <- lm(y_clr_subA ~ X); m_clr_subB <- lm(y_clr_subB ~ X)
m_plr_full <- lm(y_plr_full ~ X); m_plr_subA <- lm(y_plr_subA ~ X); m_plr_subB <- lm(y_plr_subB ~ X)


# =============================================================================
# 7. TABLE 4.2
# =============================================================================

get_slope <- function(model) coef(model)[["X"]]

table_4_2 <- data.frame(
  Method = c("Raw", "CLR", "PLR"),
  Full   = c(get_slope(m_raw_full), get_slope(m_clr_full), get_slope(m_plr_full)),
  SubA   = c(get_slope(m_raw_subA), get_slope(m_clr_subA), get_slope(m_plr_subA)),
  SubB   = c(get_slope(m_raw_subB), get_slope(m_clr_subB), get_slope(m_plr_subB))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# ---- TABLE 4.2 OUTPUT ----
print(table_4_2)
write.csv(table_4_2, "table_4_2_subcomp_results.csv", row.names = FALSE)
cat("Table 4.2 saved to table_4_2_subcomp_results.csv\n")


# =============================================================================
# 8. FIGURE 4.3 — SUBCOMPOSITIONAL COHERENCE PLOT
# =============================================================================
# Lines connecting Full, SubA, SubB show how each method's estimate changes
# (or doesn't) when the compositional context changes.
# PLR line is flat — this is subcompositional coherence.

figure_4_3_df <- data.frame(
  Method      = rep(c("Raw", "CLR", "PLR"), each = 3),
  Composition = rep(c("Full", "SubA", "SubB"), times = 3),
  Estimate    = c(
    get_slope(m_raw_full), get_slope(m_raw_subA), get_slope(m_raw_subB),
    get_slope(m_clr_full), get_slope(m_clr_subA), get_slope(m_clr_subB),
    get_slope(m_plr_full), get_slope(m_plr_subA), get_slope(m_plr_subB)
  )
) %>%
  mutate(Composition = factor(Composition, levels = c("Full", "SubA", "SubB")))

fig_4_3 <- ggplot(figure_4_3_df,
                  aes(x = Composition, y = Estimate, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = round(Estimate, 2)), vjust = -0.8, size = 4) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") +
  labs(title = "Subcompositional Coherence: Estimated Effect Across Compositional Contexts",
       x = "Compositional Context",
       y = "Estimated Effect of Binary Predictor") +
  theme_minimal(base_size = 14) +
  theme(strip.text    = element_text(face = "bold"),
        plot.title    = element_text(face = "bold", hjust = 0.5))

# ---- FIGURE 4.3 OUTPUT ----
print(fig_4_3)
ggsave("figure_4_3_subcomp_coherence.png", fig_4_3, width = 10, height = 4.5, dpi = 300)
cat("Figure 4.3 saved to figure_4_3_subcomp_coherence.png\n")
