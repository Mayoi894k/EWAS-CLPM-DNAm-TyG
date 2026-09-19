library(data.table)
library(lavaan)
library(MASS)

# ============================================================
# 0) Extract autoregressive coefficients from real CLPM results
# ============================================================
clpm_all <- fread("path/to/data/CLPM_lavaan_allCpGs_results_adj_BMI_age_T0T1_v2.csv")

message("===== Real Autoregressive Coefficient Distribution =====")
message("ar_M (beta_M0_to_M1):")
print(summary(clpm_all$beta_M0_to_M1))
message("ar_TyG (beta_TyG0_to_TyG1):")
print(summary(clpm_all$beta_TyG0_to_TyG1))

ar_M_real   <- median(clpm_all$beta_M0_to_M1,     na.rm = TRUE)
ar_TyG_real <- median(clpm_all$beta_TyG0_to_TyG1, na.rm = TRUE)
message(sprintf("Real ar_M   median: %.3f", ar_M_real))
message(sprintf("Real ar_TyG median: %.3f", ar_TyG_real))

# ============================================================
# 1) Single simulation function (unchanged core logic)
# ============================================================
simulate_once <- function(n           = 100,
                           true_effect = 0.3,
                           ar_M        = 0.6,
                           ar_TyG      = 0.6,
                           noise_sd    = 0.5) {
  Sigma  <- matrix(c(1, true_effect, true_effect, 1), 2, 2)
  T0     <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  M_T0   <- T0[, 1]
  TyG_T0 <- T0[, 2]

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

  if (is.null(fit) || !lavInspect(fit, "converged")) {
    return(data.frame(
      ewas_beta = ewas_beta,
      ewas_p    = ewas_p,
      clpm_beta = NA_real_,
      clpm_p    = NA_real_,
      converged = FALSE
    ))
  }

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

# ============================================================
# 2) Main simulation: 500 iterations using real parameters
# ============================================================
message("\n===== Part A: 500-iteration simulation with real parameters =====")
set.seed(42)
n_iter <- 500

main_results <- rbindlist(lapply(1:n_iter, function(i) {
  simulate_once(
    n           = 100,
    true_effect = 0.3,
    ar_M        = ar_M_real,
    ar_TyG      = ar_TyG_real
  )
}))

main_results <- main_results[converged == TRUE]
main_results[, sign_opposite := sign(ewas_beta) != sign(clpm_beta)]

prop_opp_main <- mean(main_results$sign_opposite)
message(sprintf("Valid iterations: %d / %d", nrow(main_results), n_iter))
message(sprintf("Proportion with opposite signs (real parameters, true_effect=0.3): %.1f%%",
                 prop_opp_main * 100))

# ============================================================
# 2b) Zero-effect scenario (true_effect = 0) under real
#     autoregressive parameters — key reviewer-requested scenario
# ============================================================
message("\n===== Part A2: Zero-effect scenario (true_effect = 0) =====")
set.seed(42)

zero_effect_results <- rbindlist(lapply(1:n_iter, function(i) {
  simulate_once(
    n           = 100,
    true_effect = 0,
    ar_M        = ar_M_real,
    ar_TyG      = ar_TyG_real
  )
}))

zero_effect_results <- zero_effect_results[converged == TRUE]
zero_effect_results[, sign_opposite := sign(ewas_beta) != sign(clpm_beta)]

prop_opp_zero <- mean(zero_effect_results$sign_opposite)
message(sprintf("Valid iterations: %d / %d", nrow(zero_effect_results), n_iter))
message(sprintf("Proportion with opposite signs (true_effect = 0, null): %.1f%%",
                 prop_opp_zero * 100))
message("If this proportion is close to 50%%, it supports the interpretation that\n",
        "sign discordance under non-null true effects reflects a structural/statistical\n",
        "artifact of comparing a difference-score EWAS coefficient to a conditional\n",
        "autoregressive CLPM path coefficient, rather than genuine reversed temporal signal.")

# ============================================================
# 3) Sensitivity analysis: parameter grid (now includes true_effect = 0)
# ============================================================
message("\n===== Part B: Sensitivity Analysis (parameter grid) =====")
set.seed(123)

grid <- expand.grid(
  ar_M        = c(0.2, 0.4, 0.6, 0.8),
  ar_TyG      = c(0.2, 0.4, 0.6, 0.8),
  true_effect = c(0, 0.1, 0.3, 0.5)   # <-- added 0 (null) scenario
)

grid_results <- rbindlist(lapply(1:nrow(grid), function(i) {
  ar_M_i   <- grid$ar_M[i]
  ar_TyG_i <- grid$ar_TyG[i]
  eff_i    <- grid$true_effect[i]

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

message("\n===== Grid results at real ar_M/ar_TyG across all true_effect settings =====")
real_setting <- grid_results[
  abs(ar_M   - round(ar_M_real,   1)) < 1e-6 &
  abs(ar_TyG - round(ar_TyG_real, 1)) < 1e-6
]
print(real_setting)

# ============================================================
# 4) Permutation / random-CpG negative control for the
#    FULL nomination-plus-CLPM procedure
# ============================================================
message("\n===== Part C: Permutation / random-CpG negative control =====")
set.seed(2024)

n_subjects   <- 100
n_cpgs       <- 1000
n_replicates <- 100
nom_p_cutoff <- 1e-4
fdr_cutoff   <- 0.05

run_negative_control_once <- function(n_subjects, n_cpgs, ar_M, ar_TyG,
                                       nom_p_cutoff, fdr_cutoff, noise_sd = 0.5) {

  TyG_T0 <- rnorm(n_subjects, 0, 1)
  TyG_T1 <- ar_TyG * TyG_T0 + rnorm(n_subjects, 0, noise_sd)
  delta_TyG <- TyG_T1 - TyG_T0

  ewas_res <- rbindlist(lapply(1:n_cpgs, function(k) {
    M_T0 <- rnorm(n_subjects, 0, 1)
    M_T1 <- ar_M * M_T0 + rnorm(n_subjects, 0, noise_sd)
    delta_M <- M_T1 - M_T0

    fit <- lm(delta_M ~ delta_TyG)
    sm  <- summary(fit)$coefficients

    data.frame(
      cpg_id    = k,
      ewas_beta = as.numeric(sm["delta_TyG", "Estimate"]),
      ewas_p    = as.numeric(sm["delta_TyG", "Pr(>|t|)"]),
      M_T0      = I(list(M_T0)),
      M_T1      = I(list(M_T1))
    )
  }))

  nominated <- ewas_res[ewas_p < nom_p_cutoff]
  n_nominated <- nrow(nominated)

  if (n_nominated == 0) {
    return(data.frame(n_nominated = 0, n_clpm_fdr_sig = 0, prop_false_positive = NA_real_))
  }

  clpm_res <- rbindlist(lapply(seq_len(n_nominated), function(idx) {
    M_T0 <- nominated$M_T0[[idx]]
    M_T1 <- nominated$M_T1[[idx]]
    d <- data.frame(M_T0, M_T1, TyG_T0, TyG_T1)

    model <- '
      M_T1   ~ a*M_T0   + c*TyG_T0
      TyG_T1 ~ b*TyG_T0 + d*M_T0
      M_T0   ~~ TyG_T0
      M_T1   ~~ TyG_T1
    '
    fit <- tryCatch(sem(model, data = d, estimator = "ML"), error = function(e) NULL)
    if (is.null(fit) || !lavInspect(fit, "converged")) {
      return(data.frame(clpm_p = NA_real_))
    }
    pe <- parameterEstimates(fit)
    clpm_row <- pe[pe$op == "~" & pe$lhs == "TyG_T1" & pe$rhs == "M_T0", ]
    data.frame(clpm_p = as.numeric(clpm_row$pvalue[1]))
  }))

  clpm_res <- clpm_res[!is.na(clpm_p)]
  clpm_res[, clpm_fdr := p.adjust(clpm_p, method = "BH")]
  n_clpm_fdr_sig <- sum(clpm_res$clpm_fdr < fdr_cutoff, na.rm = TRUE)

  data.frame(
    n_nominated         = n_nominated,
    n_clpm_fdr_sig       = n_clpm_fdr_sig,
    prop_false_positive  = n_clpm_fdr_sig / n_nominated
  )
}

negctrl_results <- rbindlist(lapply(1:n_replicates, function(r) {
  run_negative_control_once(
    n_subjects   = n_subjects,
    n_cpgs       = n_cpgs,
    ar_M         = ar_M_real,
    ar_TyG       = ar_TyG_real,
    nom_p_cutoff = nom_p_cutoff,
    fdr_cutoff   = fdr_cutoff
  )
}))

message(sprintf(
  "Negative control (null CpGs, full nomination+CLPM pipeline), %d replicates:",
  n_replicates
))
message(sprintf("Mean number of CpGs nominated per replicate (P<%.0e): %.2f",
                 nom_p_cutoff, mean(negctrl_results$n_nominated)))
message(sprintf("Mean number of CLPM FDR<%.2f 'significant' false positives per replicate: %.3f",
                 fdr_cutoff, mean(negctrl_results$n_clpm_fdr_sig)))
message(sprintf("Mean false-positive proportion among nominated CpGs: %.3f (SD %.3f)",
                 mean(negctrl_results$prop_false_positive, na.rm = TRUE),
                 sd(negctrl_results$prop_false_positive, na.rm = TRUE)))

# ============================================================
# 5) Save all results
# ============================================================
out_dir <- "path/to/output/directory"   # <-- output directory placeholder
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(main_results,        file.path(out_dir, "sim_main_results.csv"))
fwrite(zero_effect_results, file.path(out_dir, "sim_zero_effect_results.csv"))
fwrite(grid_results,        file.path(out_dir, "sim_sensitivity_grid.csv"))
fwrite(negctrl_results,     file.path(out_dir, "sim_negative_control_results.csv"))

message("\nSaved: sim_main_results.csv")
message("Saved: sim_zero_effect_results.csv")
message("Saved: sim_sensitivity_grid.csv (now includes true_effect = 0)")
message("Saved: sim_negative_control_results.csv")
message("Output directory: ", out_dir)

# ============================================================
# 6) Figures
# ============================================================
out_fig <- file.path(out_dir, "simulation_plots.pdf")
pdf(out_fig, width = 12, height = 8)

par(mfrow = c(2, 3))

hist(main_results$ewas_beta,
     main = sprintf("EWAS beta (ar_M=%.2f, ar_TyG=%.2f, true_effect=0.3)", ar_M_real, ar_TyG_real),
     xlab = "beta", col = "steelblue", border = "white", breaks = 30)
abline(v = 0, col = "red", lwd = 2, lty = 2)

hist(main_results$clpm_beta,
     main = "CLPM beta (M_T0 -> TyG_T1, true_effect=0.3)",
     xlab = "beta", col = "darkorange", border = "white", breaks = 30)
abline(v = 0, col = "red", lwd = 2, lty = 2)

plot(main_results$ewas_beta, main_results$clpm_beta,
     xlab = "EWAS beta", ylab = "CLPM beta",
     main = sprintf("EWAS vs CLPM (opposite signs: %.1f%%)", prop_opp_main * 100),
     pch = 16, col = rgb(0.2, 0.4, 0.8, 0.4), cex = 0.8)
abline(h = 0, col = "gray", lty = 2)
abline(v = 0, col = "gray", lty = 2)

plot(zero_effect_results$ewas_beta, zero_effect_results$clpm_beta,
     xlab = "EWAS beta", ylab = "CLPM beta",
     main = sprintf("Null scenario (true_effect=0): opposite signs %.1f%%", prop_opp_zero * 100),
     pch = 16, col = rgb(0.6, 0.6, 0.6, 0.4), cex = 0.8)
abline(h = 0, col = "gray", lty = 2)
abline(v = 0, col = "gray", lty = 2)

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
      main = "Proportion opposite signs (true_effect = 0.3)",
      col  = heat.colors(20))
points(ar_M_real, ar_TyG_real, pch = 4, cex = 2, lwd = 3, col = "blue")
legend("topleft", legend = "Real parameters", pch = 4, col = "blue", bty = "n")

grid_sub0 <- grid_results[true_effect == 0]
mat0 <- matrix(NA, 4, 4)
for (r in 1:nrow(grid_sub0)) {
  i <- which(ar_vals == grid_sub0$ar_M[r])
  j <- which(ar_vals == grid_sub0$ar_TyG[r])
  mat0[i, j] <- grid_sub0$prop_opposite[r]
}
image(ar_vals, ar_vals, mat0,
      xlab = "ar_M", ylab = "ar_TyG",
      main = "Proportion opposite signs (true_effect = 0, null)",
      col  = heat.colors(20))
points(ar_M_real, ar_TyG_real, pch = 4, cex = 2, lwd = 3, col = "blue")
legend("topleft", legend = "Real parameters", pch = 4, col = "blue", bty = "n")

dev.off()
message("Saved figure: ", out_fig)

# ============================================================
# 7) Summary
# ============================================================
message("\n========== Simulation Summary ==========")
message(sprintf("Real autoregressive coefficients: ar_M=%.3f, ar_TyG=%.3f", ar_M_real, ar_TyG_real))
message(sprintf("Proportion opposite signs (true_effect=0.3, real params): %.1f%%", prop_opp_main * 100))
message(sprintf("Proportion opposite signs (true_effect=0, null, real params): %.1f%%", prop_opp_zero * 100))
message("Sensitivity analysis across true_effect settings at real ar_M/ar_TyG:")
print(real_setting[, .(ar_M, ar_TyG, true_effect, prop_opposite)])
message("\nNegative control (full nomination + CLPM pipeline under global null):")
message(sprintf("Mean CpGs nominated per replicate: %.2f", mean(negctrl_results$n_nominated)))
message(sprintf("Mean CLPM FDR<0.05 false positives per replicate: %.3f",
                 mean(negctrl_results$n_clpm_fdr_sig)))
message("=========================================")