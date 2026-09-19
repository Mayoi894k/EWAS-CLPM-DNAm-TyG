library(data.table)

# ============================================================
# Compute three diagnostic correlation coefficients across all CpG sites:
#   (1) cor(M_T0, TyG_T0)  -- baseline cross-sectional correlation
#   (2) cor(delta_M, delta_TyG) -- EWAS estimand (difference-score correlation)
#   (3) cor(M_T0, TyG_T1)  -- CLPM estimand (prospective cross-lagged correlation)
# ============================================================

# ---------------------------
# 1. Path settings
# ---------------------------
in_file  <- "path/to/data/CLPM_input.csv"
out_dir  <- "path/to/output/directory"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------------------------
# 2. Read data
# ---------------------------
dt <- fread(in_file)

cat("Number of rows (subjects):", nrow(dt), "\n")
cat("Number of columns:", ncol(dt), "\n")

# Check whether TyG-related columns exist
req_cols <- c("T0_TyG", "T1_TyG")
missing_req <- setdiff(req_cols, names(dt))
if (length(missing_req) > 0) {
  stop("Required columns are missing from the input file: ", paste(missing_req, collapse = ", "))
}

dt[, T0_TyG := as.numeric(T0_TyG)]
dt[, T1_TyG := as.numeric(T1_TyG)]
dt[, delta_TyG := T1_TyG - T0_TyG]

# ---------------------------
# 3. Detect all CpG sites (paired by _T0 / _T1 suffix)
# ---------------------------
cpg_t0_cols <- grep("_T0$", names(dt), value = TRUE)
cpg_t1_cols <- grep("_T1$", names(dt), value = TRUE)

# Exclude TyG itself (avoid mistakenly identifying it as a CpG)
cpg_t0_cols <- setdiff(cpg_t0_cols, "T0_TyG")
cpg_t1_cols <- setdiff(cpg_t1_cols, "T1_TyG")

cpg_base_t0 <- sub("_T0$", "", cpg_t0_cols)
cpg_base_t1 <- sub("_T1$", "", cpg_t1_cols)
common_cpg  <- intersect(cpg_base_t0, cpg_base_t1)

cat("Number of CpG sites detected:", length(common_cpg), "\n")

if (length(common_cpg) == 0) {
  stop("No CpG sites were detected. Please check whether column name suffixes are _T0 / _T1.")
}

# ---------------------------
# 4. Compute the three correlation coefficients for each CpG site
# ---------------------------
compute_cor_one_cpg <- function(cpg) {
  col0 <- paste0(cpg, "_T0")
  col1 <- paste0(cpg, "_T1")

  M_T0 <- as.numeric(dt[[col0]])
  M_T1 <- as.numeric(dt[[col1]])
  delta_M <- M_T1 - M_T0

  # Retain complete observations required for each correlation (pairwise complete)
  cor_MT0_TyGT0 <- suppressWarnings(cor(M_T0, dt$T0_TyG, use = "pairwise.complete.obs"))
  cor_dM_dTyG   <- suppressWarnings(cor(delta_M, dt$delta_TyG, use = "pairwise.complete.obs"))
  cor_MT0_TyGT1 <- suppressWarnings(cor(M_T0, dt$T1_TyG, use = "pairwise.complete.obs"))

  # Effective sample sizes (using the pairwise-complete N for each pair)
  n_MT0_TyGT0 <- sum(complete.cases(M_T0, dt$T0_TyG))
  n_dM_dTyG   <- sum(complete.cases(delta_M, dt$delta_TyG))
  n_MT0_TyGT1 <- sum(complete.cases(M_T0, dt$T1_TyG))

  data.frame(
    unit               = cpg,
    cor_MT0_TyGT0      = cor_MT0_TyGT0,
    n_MT0_TyGT0        = n_MT0_TyGT0,
    cor_dM_dTyG        = cor_dM_dTyG,
    n_dM_dTyG          = n_dM_dTyG,
    cor_MT0_TyGT1      = cor_MT0_TyGT1,
    n_MT0_TyGT1        = n_MT0_TyGT1,
    stringsAsFactors   = FALSE
  )
}

results <- rbindlist(lapply(common_cpg, compute_cor_one_cpg))

# ---------------------------
# 5. Flag whether signs agree/disagree (EWAS estimand vs. CLPM estimand direction)
# ---------------------------
results[, sign_opposite_dMdTyG_vs_MT0TyGT1 := sign(cor_dM_dTyG) != sign(cor_MT0_TyGT1)]

cat("\n===== Summary =====\n")
cat("Total number of CpGs:", nrow(results), "\n")
cat("Number of sites where cor(delta_M, delta_TyG) and cor(M_T0, TyG_T1) have opposite signs:",
    sum(results$sign_opposite_dMdTyG_vs_MT0TyGT1, na.rm = TRUE), "\n")
cat("Proportion with opposite signs:",
    round(mean(results$sign_opposite_dMdTyG_vs_MT0TyGT1, na.rm = TRUE) * 100, 1), "%\n")

# ---------------------------
# 6. Save results
# ---------------------------
out_path <- file.path(out_dir, "diagnostic_correlations_allCpGs.csv")
fwrite(results, out_path)
cat("\nResults saved to:", out_path, "\n")

# Also save a version sorted by the absolute value of cor_dM_dTyG for convenience
results_sorted <- results[order(-abs(cor_dM_dTyG))]
out_path_sorted <- file.path(out_dir, "diagnostic_correlations_allCpGs_sorted.csv")
fwrite(results_sorted, out_path_sorted)
cat("Sorted version saved to:", out_path_sorted, "\n")