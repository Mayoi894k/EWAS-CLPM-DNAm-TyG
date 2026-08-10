options(stringsAsFactors = FALSE)

library(data.table)
library(lme4)
library(lmerTest)
library(parallel)

# -------------------------
# Config
# -------------------------
# Path to the merged wide-format input CSV file (longitudinal data)
in_file <- "/path/to/your/input_data.csv"

# Output directory for all result files
out_dir <- "/path/to/your/output_directory"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum number of complete observations required per CpG site to fit the model
min_obs_per_site <- 6

# Use up to 18 cores (leaves 1 core free for the OS)
n_cores <- min(detectCores() - 1, 18)
if (n_cores < 1) n_cores <- 1

# P-value thresholds for saving stratified result files
p_cutoffs <- c(0.0001, 0.001, 0.005, 0.01, 0.05, 1)

# Benjamini-Hochberg FDR threshold for the significant-hits output file
fdr_cutoff <- 0.05

# -------------------------
# 1) Read merged wide data
# -------------------------
# Expected format: one row per subject; baseline columns suffixed with .x,
# follow-up columns suffixed with .y; CpG methylation columns suffixed with _T0/_T1.
message("Reading: ", in_file)
dt <- fread(in_file)

# ------------------------------------------------------------------
# Required non-CpG columns that must exist in the input file:
#   subject_id   - unique subject identifier
#   family_ID.x  - family ID at baseline (family_ID.y used as fallback)
#   T0_TyG       - baseline TyG index
#   T1_TyG       - follow-up TyG index
#   T0_age       - baseline age
#   T0_gender    - baseline sex (1 = male, 2 = female)
# ------------------------------------------------------------------
req_cols <- c("subject_id", "family_ID.x", "T0_TyG", "T1_TyG", "T0_age", "T0_gender")
missing_req <- setdiff(req_cols, names(dt))
if (length(missing_req) > 0) {
  stop("Missing required columns in input file: ", paste(missing_req, collapse = ", "))
}

# ------------------------------------------------------------------
# Build family_ID: use baseline value; fall back to follow-up if missing
# ------------------------------------------------------------------
family_vec <- dt[["family_ID.x"]]
if ("family_ID.y" %in% names(dt)) {
  miss <- is.na(family_vec) | family_vec == ""
  if (any(miss)) family_vec[miss] <- dt[["family_ID.y"]][miss]
}
dt[, family_ID := factor(as.character(family_vec))]

# Coerce covariates and exposure to numeric
dt[, T0_TyG      := as.numeric(T0_TyG)]
dt[, T1_TyG      := as.numeric(T1_TyG)]
dt[, T0_age      := as.numeric(T0_age)]
dt[, T0_gender   := as.numeric(T0_gender)]

# Recode sex: 1 (male) -> 0, 2 (female) -> 1; other values -> NA
dt[, T0_gender01 := fifelse(T0_gender == 1, 0,
                     fifelse(T0_gender == 2, 1, NA_real_))]

# Compute change in TyG index (exposure variable for the EWAS)
dt[, delta_TyG := T1_TyG - T0_TyG]

# -------------------------
# 2) Detect CpG columns and build CpG pairs
# -------------------------
# CpG methylation columns are expected to be named <CpGid>_T0 and <CpGid>_T1
cpg_t0_cols <- grep("_T0$", names(dt), value = TRUE)
cpg_t1_cols <- grep("_T1$", names(dt), value = TRUE)

if (length(cpg_t0_cols) == 0 || length(cpg_t1_cols) == 0) {
  stop("No CpG columns found with suffix _T0/_T1. Please verify input column naming.")
}

cpg_base_t0 <- sub("_T0$", "", cpg_t0_cols)
cpg_base_t1 <- sub("_T1$", "", cpg_t1_cols)
common_cpg  <- intersect(cpg_base_t0, cpg_base_t1)
if (length(common_cpg) == 0) stop("No matched CpG pairs between _T0 and _T1 columns.")
message("Detected CpG pairs: ", length(common_cpg))

# -------------------------
# 3) Per-CpG delta EWAS model
#
#    Model: delta_M ~ delta_TyG + T0_age + T0_gender01 + (1 | family_ID)
#      delta_M     = change in methylation (T1 - T0) at each CpG
#      delta_TyG   = change in TyG index (primary exposure)
#      T0_age      = baseline age (covariate)
#      T0_gender01 = baseline sex, 0/1 coded (covariate)
#      (1|family_ID) = random intercept for family clustering
# -------------------------
fit_one_cpg <- function(i) {
  cpg  <- common_cpg[i]
  col0 <- paste0(cpg, "_T0")
  col1 <- paste0(cpg, "_T1")

  # Compute per-CpG change in M-value (or beta-value, depending on input)
  delta_M <- as.numeric(dt[[col1]]) - as.numeric(dt[[col0]])

  df <- data.frame(
    unit        = cpg,
    delta_M     = delta_M,
    delta_TyG   = dt$delta_TyG,
    T0_age      = dt$T0_age,
    T0_gender01 = dt$T0_gender01,
    family_ID   = dt$family_ID,
    stringsAsFactors = FALSE
  )

  # Keep only rows with complete data for all model variables
  df <- df[
    is.finite(df$delta_M) &
      is.finite(df$delta_TyG) &
      is.finite(df$T0_age) &
      is.finite(df$T0_gender01) &
      !is.na(df$family_ID),
    ,
    drop = FALSE
  ]

  Nobs <- nrow(df)

  # Skip CpG if too few observations or only one family (random effect not estimable)
  if (Nobs < min_obs_per_site || length(unique(df$family_ID)) < 2) {
    return(data.frame(
      unit    = cpg,
      Estimate = NA_real_, StdErr  = NA_real_,
      t_value  = NA_real_, P.Value = NA_real_,
      Nobs     = Nobs,
      stringsAsFactors = FALSE
    ))
  }

  out <- tryCatch({
    fit <- lmer(
      delta_M ~ delta_TyG + T0_age + T0_gender01 + (1 | family_ID),
      data = df, REML = FALSE
    )
    sm <- summary(fit)$coefficients

    if (!"delta_TyG" %in% rownames(sm)) {
      return(data.frame(
        unit    = cpg,
        Estimate = NA_real_, StdErr  = NA_real_,
        t_value  = NA_real_, P.Value = NA_real_,
        Nobs     = Nobs
      ))
    }

    est  <- as.numeric(sm["delta_TyG", "Estimate"])
    se   <- as.numeric(sm["delta_TyG", "Std. Error"])
    tval <- as.numeric(sm["delta_TyG", "t value"])
    pcol <- if ("Pr(>|t|)" %in% colnames(sm)) "Pr(>|t|)" else NA_character_
    pval <- if (!is.na(pcol)) as.numeric(sm["delta_TyG", pcol]) else NA_real_

    data.frame(
      unit     = cpg,
      Estimate = est,
      StdErr   = se,
      t_value  = tval,
      P.Value  = pval,
      Nobs     = Nobs,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    # Return NA row on convergence or other model failures
    data.frame(
      unit    = cpg,
      Estimate = NA_real_, StdErr  = NA_real_,
      t_value  = NA_real_, P.Value = NA_real_,
      Nobs     = Nobs,
      stringsAsFactors = FALSE
    )
  })

  out
}

message(
  "Starting adjusted delta-EWAS (lmer) over ", length(common_cpg),
  " CpGs using ", n_cores, " cores."
)
res_list <- mclapply(seq_along(common_cpg), fit_one_cpg, mc.cores = n_cores)
res_df   <- rbindlist(res_list)

# Apply Benjamini-Hochberg FDR correction across all CpGs
res_df[, adj.P.Val := p.adjust(P.Value, method = "BH")]

# -------------------------
# 4) Genomic inflation factor (lambda)
#    lambda > 1 suggests systematic inflation; lambda ~ 1 is expected under the null.
# -------------------------
pvals <- as.numeric(res_df$P.Value)
pvals <- pvals[!is.na(pvals) & is.finite(pvals) & pvals > 0 & pvals <= 1]

lambda <- NA_real_
if (length(pvals) >= 10) {
  chisq  <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
  message("\n========== EWAS Genomic Inflation (lambda) ==========")
  message("Number of p-values used: ", length(pvals))
  message(sprintf("Lambda (GC): %.4f", lambda))
  message("=====================================================\n")
} else {
  message("WARNING: Too few valid p-values to compute lambda (n=", length(pvals), ").")
}

# Save lambda summary to a plain-text file
writeLines(
  c(
    paste0("Run time: ",       format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Input: ",          in_file),
    paste0("CpGs: ",           length(common_cpg)),
    paste0("Valid p-values: ", length(pvals)),
    paste0("Lambda(GC): ",     ifelse(is.na(lambda), "NA", sprintf("%.6f", lambda)))
  ),
  con = file.path(out_dir, "lambda.txt")
)

# -------------------------
# 5) Save full results
# -------------------------
full_csv <- file.path(out_dir, "delta_lmer_full_results_adj_age_sex.csv")
fwrite(res_df, full_csv)
message("Saved full results: ", full_csv)

# -------------------------
# 6) Save FDR-significant results
# -------------------------
sig_fdr     <- res_df[!is.na(adj.P.Val) & adj.P.Val < fdr_cutoff]
sig_fdr     <- sig_fdr[order(P.Value)]
sig_fdr_csv <- file.path(out_dir, paste0("sig_cpgs_FDR", fdr_cutoff, "_adj_age_sex.csv"))
fwrite(sig_fdr, sig_fdr_csv)
message("Saved FDR<", fdr_cutoff, ": ", nrow(sig_fdr), " -> ", sig_fdr_csv)

# -------------------------
# 7) Save results stratified by nominal p-value thresholds
# -------------------------
for (pcut in p_cutoffs) {
  if (pcut >= 1) {
    sig     <- res_df
    out_csv <- file.path(out_dir, "sig_cpgs_p1_adj_age_sex.csv")
  } else {
    sig <- res_df[!is.na(P.Value) & P.Value < pcut]
    if (abs(pcut - 0.0001) < 1e-12) {
      out_csv <- file.path(out_dir, "sig_cpgs_p1e-04_adj_age_sex.csv")
    } else {
      out_csv <- file.path(out_dir, paste0("sig_cpgs_p", pcut, "_adj_age_sex.csv"))
    }
  }
  sig <- sig[order(P.Value)]
  fwrite(sig, out_csv)
  message("Saved p<", pcut, ": ", nrow(sig), " -> ", out_csv)
}

message("\n==========================================")
message("Adjusted delta-EWAS (lmer) finished successfully!")
message("Output directory: ", out_dir)
message("==========================================")