library(data.table)
library(lavaan)
library(MASS)

# -------------------------
# 0) Extract autoregressive coefficients from real CLPM results
# -------------------------
clpm_all <- fread("path/to/CLPM_lavaan_allCpGs_results_adj_BMI_age_T0T1.csv")

# Inspect the distribution of real autoregressive coefficients
message("===== Real Autoregressive Coefficient Distribution =====")
message("ar_M (beta_M0_to_M1):")
print(summary(clpm_all$beta_M0_to_M1))
message("ar_TyG (beta_TyG0_to_TyG1):")
print(summary(clpm_all$beta_TyG0_to_TyG1))

# Use medians as reference parameters for simulation
ar_M_real   <- median(clpm_all$beta_M0_to_M1,     na.rm = TRUE)
ar_TyG_real <- median(clpm_all$beta_TyG0_to_TyG1, na.rm = TRUE)
message(sprintf("Real ar_M   median: %.3f", ar_M_real))
message(sprintf("Real ar_TyG median: %.3f", ar_TyG_real))

# -------------------------
# 1) Single simulation function
# -------------------------
simulate_once <- function(n           = 100,
                          true_effect = 0.3,
                          ar_M        = 0.6,
                          ar_TyG      = 0.6,
                          noise_sd    = 0.5) {
  # Generate T0 joint distribution with specified correlation
  Sigma  <- matrix(c(1, true_effect, true_effect, 1), 2, 2)
  T0     <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  M_T0   <- T0[, 1]
  TyG_T0 <- T0[, 2]

  # Generate T1 variables using cross-lagged structure
  M_T1   <- ar_M   * M_T0   + true_effect * TyG_T0 + rnorm(n, 0, noise_sd)
  TyG_T1 <- ar_TyG * TyG_T0 + true_effect * M_T0   + rnorm(n, 0, noise_sd)

  d <- data.frame(M_T0, M_T1, TyG_T0, TyG_T1)

  # ---------- EWAS: regress change in M on change in TyG ----------
  delta_M   <- M_T1 - M_T0
  delta_TyG <- TyG_T1 - TyG_T0
  ewas_fit  <- lm(delta_M ~ delta_TyG)
  ewas_beta <- as.numeric(coef(ewas_fit)["delta_TyG"])
  ewas_p    <- as.numeric(summary(ewas_fit)$coefficients["delta_TyG", "Pr(>|t|)"])

  # ---------- CLPM: cross-lagged panel model via lavaan ----------
  model <- '
    M_T1   ~ a*M_T0   + c*TyG_T0
    TyG_T1 ~ b*TyG_T0 + d*M_T0
    M_T0   ~~ TyG_T0
    M_T1   ~~ TyG_T1
  '
  fit <- tryCatch(
    sem(model, data = d, estimator = "ML"),
    error = function(e) NULL
  )

  # Return EWAS results only if CLPM did not converge
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    return(data.frame(
      ewas_beta = ewas_beta,
      ewas_p    = ewas_p,
      clpm_beta = NA_real_,
      clpm_p    = NA_real_,
      converged = FALSE
    ))
  }

  # Extract the cross-lagged path estimate: M_T0 -> TyG_T1
  pe        <- parameterEstimates(fit)
  clpm_row  <- pe[pe$op == "~" & pe$lhs == "TyG_T1" & pe$rhs == "M_T0", ]

  data.frame(
    ewas_beta = ewas_beta,
    ewas_p    = ewas_p,
    clpm_beta = as.numeric(clpm_row$est[1]),
    clpm_p    = as.numeric(clpm_row$pvalue[1]),
    converged = TRUE
  )
}

# -------------------------
# 2) Main simulation: 500 iterations using real parameters
# -------------------------
message("\n===== Part A: 500-iteration simulation with real parameters =====")
set.seed(42)
n_iter <- 500

main_results <- rbindlist(lapply(1:n_iter, function(i) {
  simulate_once(
    n           = 100,
    true_effect = 0.3,       # moderate effect size
    ar_M        = ar_M_real,
    ar_TyG      = ar_TyG_real
  )
}))

# Keep only converged runs and flag sign discordance between EWAS and CLPM
main_results <- main_results[converged == TRUE]
main_results[, sign_opposite := sign(ewas_beta) != sign(clpm_beta)]

prop_opp_main <- mean(main_results$sign_opposite)
message(sprintf("Valid iterations: %d / %d", nrow(main_results), n_iter))
message(sprintf("Proportion with opposite signs (real parameters): %.1f%%", prop_opp_main * 100))

# -------------------------
# 3) Sensitivity analysis: parameter grid
# -------------------------
message("\n===== Part B: Sensitivity Analysis (parameter grid) =====")
set.seed(123)

grid <- expand.grid(
  ar_M        = c(0.2, 0.4, 0.6, 0.8),
  ar_TyG      = c(0.2, 0.4, 0.6, 0.8),
  true_effect = c(0.1, 0.3, 0.5)
)

grid_results <- rbindlist(lapply(1:nrow(grid), function(i) {
  ar_M_i   <- grid$ar_M[i]
  ar_TyG_i <- grid$ar_TyG[i]
  eff_i    <- grid$true_effect[i]

  # Run 200 simulations per parameter combination
  sims <- rbindlist(lapply(1:200, function(j) {
    simulate_once(100, eff_i, ar_M_i, ar_TyG_i)
  }))

  sims <- sims[converged == TRUE]
  if (nrow(sims) == 0) return(NULL)

  data.frame(
    ar_M          = ar_M_i,
    ar_TyG        = ar_TyG_i,
    true_effect   = eff_i,
    n_converged   = nrow(sims),
    prop_opposite = round(mean(sign(sims$ewas_beta) != sign(sims$clpm_beta)), 3)
  )
}))

print(grid_results)

# -------------------------
# 4) Save results
# -------------------------
out_dir <- "output"   # set to your preferred output directory

fwrite(main_results,  file.path(out_dir, "sim_main_results.csv"))
fwrite(grid_results,  file.path(out_dir, "sim_sensitivity_grid.csv"))
message("\nSaved: sim_main_results.csv")
message("Saved: sim_sensitivity_grid.csv")

# -------------------------
# 5) Figures
# -------------------------
out_fig <- file.path(out_dir, "simulation_plots.pdf")
pdf(out_fig, width = 12, height = 5)

# Panel 1 & 2: Beta distributions for EWAS and CLPM
par(mfrow = c(1, 3))

hist(main_results$ewas_beta,
     main = sprintf("EWAS beta (ar_M=%.2f, ar_TyG=%.2f)", ar_M_real, ar_TyG_real),
     xlab = "beta", col = "steelblue", border = "white", breaks = 30)
abline(v = 0, col = "red", lwd = 2, lty = 2)

hist(main_results$clpm_beta,
     main = "CLPM beta (M_T0 -> TyG_T1)",
     xlab = "beta", col = "darkorange", border = "white", breaks = 30)
abline(v = 0, col = "red", lwd = 2, lty = 2)

# Panel 3: Scatter plot of EWAS vs CLPM estimates
plot(main_results$ewas_beta, main_results$clpm_beta,
     xlab = "EWAS beta", ylab = "CLPM beta",
     main = sprintf("EWAS vs CLPM (opposite signs: %.1f%%)", prop_opp_main * 100),
     pch = 16, col = rgb(0.2, 0.4, 0.8, 0.4), cex = 0.8)
abline(h = 0, col = "gray", lty = 2)
abline(v = 0, col = "gray", lty = 2)

# Panel 4: Sensitivity heatmap for true_effect = 0.3
grid_sub <- grid_results[true_effect == 0.3]
mat <- matrix(NA, 4, 4)
ar_vals <- c(0.2, 0.4, 0.6, 0.8)
for (r in 1:nrow(grid_sub)) {
  i <- which(ar_vals == grid_sub$ar_M[r])
  j <- which(ar_vals == grid_sub$ar_TyG[r])
  mat[i, j] <- grid_sub$prop_opposite[r]
}
image(ar_vals, ar_vals, mat,
      xlab = "ar_M", ylab = "ar_TyG",
      main = "Proportion with opposite signs (true_effect = 0.3)",
      col  = heat.colors(20))
# Mark the real parameter location
points(ar_M_real, ar_TyG_real, pch = 4, cex = 2, lwd = 3, col = "blue")
legend("topleft", legend = "Real parameters", pch = 4, col = "blue", bty = "n")

dev.off()
message("Saved figure: ", out_fig)

# -------------------------
# 6) Summary
# -------------------------
message("\n========== Simulation Summary ==========")
message(sprintf("Real autoregressive coefficients: ar_M=%.3f, ar_TyG=%.3f", ar_M_real, ar_TyG_real))
message(sprintf("Proportion with opposite signs (real parameters): %.1f%%", prop_opp_main * 100))
message("Sensitivity analysis (true_effect = 0.3):")
print(grid_results[true_effect == 0.3, .(ar_M, ar_TyG, prop_opposite)])
message("=========================================")