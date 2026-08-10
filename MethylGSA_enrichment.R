# ============================================================
# methylGSA Enrichment Analysis - GO / KEGG / Reactome
# Input: sig_cpgs_p1_adj_age_sex.csv
# ============================================================

# ---- 1. Install required packages --------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Upgrade BiocManager to match R 4.6
BiocManager::install(version = "3.23")

# Install methylGSA and dependencies
BiocManager::install("methylGSA", ask = FALSE)
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE)
BiocManager::install("org.Hs.eg.db", ask = FALSE)
BiocManager::install("reactome.db",  ask = FALSE)

# ---- 2. Load packages ---------------------------------------
library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ---- 3. Read data -------------------------------------------
input_file <- "sig_cpgs_p1_adj_age_sex.csv"
output_dir <- "./"

dat <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Check data structure (confirm column names)
head(dat)
cat("Column names:", colnames(dat), "\n")
cat("Total rows:", nrow(dat), "\n")

# ---- 4. Extract CpG information -----------------------------
# Column "unit" contains CpG IDs; "adj.P.Val" contains adjusted p-values

# All tested CpGs (used as background)
all_cpg <- dat$unit

# Significant CpGs (adj.P.Val < 0.05)
sig_cpg <- dat$unit[dat$adj.P.Val < 0.05]

cat("Total CpGs:", length(all_cpg), "\n")
cat("Significant CpGs (adj.P < 0.05):", length(sig_cpg), "\n")

# If too few significant CpGs, fall back to nominal p-value threshold
if (length(sig_cpg) < 10) {
  cat("Too few significant CpGs; switching to nominal P.Value < 0.05\n")
  sig_cpg <- dat$unit[dat$P.Value < 0.05]
  cat("Significant CpGs after re-filtering:", length(sig_cpg), "\n")
}

# ---- 5. Build cpg.pval vector (recommended input for methylGSA) ----
# methylglm / methylRRA require a named p-value vector (all CpGs)
cpg_pval <- setNames(dat$adj.P.Val, dat$unit)

# ---- 6. GO Enrichment Analysis ------------------------------
cat("\n========== GO Enrichment Analysis ==========\n")

## Method 1: methylglm (logistic regression, recommended)
res_GO <- tryCatch({
  methylglm(
    cpg.pval    = cpg_pval,
    array.type  = "450K",       # Change to "EPIC" if using EPIC array
    group       = "all",        # Use all CpGs as background
    GS.type     = "GO",
    minsize     = 5,
    maxsize     = 500
  )
}, error = function(e) {
  cat("methylglm GO failed; switching to methylgometh\n", conditionMessage(e), "\n")
  methylgometh(
    sig.cpg    = sig_cpg,
    all.cpg    = all_cpg,
    array.type = "450K",
    group      = "all",
    GS.type    = "GO",
    minsize    = 5,
    maxsize    = 500
  )
})

cat("GO result rows:", nrow(res_GO), "\n")
print(head(res_GO[order(res_GO$pvalue), ], 10))

# Save GO results
go_out <- file.path(output_dir, "GO_enrichment_results.csv")
write.csv(res_GO[order(res_GO$pvalue), ], go_out, row.names = FALSE)
cat("GO results saved to:", go_out, "\n")

# ---- 7. KEGG Enrichment Analysis ----------------------------
cat("\n========== KEGG Enrichment Analysis ==========\n")

res_KEGG <- tryCatch({
  methylglm(
    cpg.pval    = cpg_pval,
    array.type  = "450K",
    group       = "all",
    GS.type     = "KEGG",
    minsize     = 5,
    maxsize     = 500
  )
}, error = function(e) {
  cat("methylglm KEGG failed; switching to methylgometh\n", conditionMessage(e), "\n")
  methylgometh(
    sig.cpg    = sig_cpg,
    all.cpg    = all_cpg,
    array.type = "450K",
    group      = "all",
    GS.type    = "KEGG",
    minsize    = 5,
    maxsize    = 500
  )
})

cat("KEGG result rows:", nrow(res_KEGG), "\n")
print(head(res_KEGG[order(res_KEGG$pvalue), ], 10))

kegg_out <- file.path(output_dir, "KEGG_enrichment_results.csv")
write.csv(res_KEGG[order(res_KEGG$pvalue), ], kegg_out, row.names = FALSE)
cat("KEGG results saved to:", kegg_out, "\n")

# ---- 8. Reactome Enrichment Analysis ------------------------
cat("\n========== Reactome Enrichment Analysis ==========\n")

res_Reactome <- tryCatch({
  methylglm(
    cpg.pval    = cpg_pval,
    array.type  = "450K",
    group       = "all",
    GS.type     = "Reactome",
    minsize     = 5,
    maxsize     = 500
  )
}, error = function(e) {
  cat("methylglm Reactome failed; switching to methylgometh\n", conditionMessage(e), "\n")
  methylgometh(
    sig.cpg    = sig_cpg,
    all.cpg    = all_cpg,
    array.type = "450K",
    group      = "all",
    GS.type    = "Reactome",
    minsize    = 5,
    maxsize    = 500
  )
})

cat("Reactome result rows:", nrow(res_Reactome), "\n")
print(head(res_Reactome[order(res_Reactome$pvalue), ], 10))

reactome_out <- file.path(output_dir, "Reactome_enrichment_results.csv")
write.csv(res_Reactome[order(res_Reactome$pvalue), ], reactome_out, row.names = FALSE)
cat("Reactome results saved to:", reactome_out, "\n")

# ---- 9. Summarize significant enrichment results (padj < 0.05) ----
cat("\n========== Summary of Significant Results (padj < 0.05) ==========\n")

summarize_sig <- function(res, name) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  # methylGSA result column may be named "padj", "BH", or "p.adjust"
  padj_col <- intersect(c("padj", "BH", "p.adjust"), colnames(res))
  if (length(padj_col) == 0) {
    cat(name, ": adjusted p-value column not found; showing results with p-value < 0.05\n")
    sig <- res[res$pvalue < 0.05, ]
  } else {
    sig <- res[res[[padj_col[1]]] < 0.05, ]
  }
  cat(name, "significant terms:", nrow(sig), "\n")
  sig$database <- name
  return(sig)
}

sig_GO       <- summarize_sig(res_GO,       "GO")
sig_KEGG     <- summarize_sig(res_KEGG,     "KEGG")
sig_Reactome <- summarize_sig(res_Reactome, "Reactome")

# Merge and save
all_sig <- do.call(rbind, Filter(Negate(is.null), list(sig_GO, sig_KEGG, sig_Reactome)))
if (!is.null(all_sig) && nrow(all_sig) > 0) {
  summary_out <- file.path(output_dir, "All_significant_enrichment_summary.csv")
  write.csv(all_sig, summary_out, row.names = FALSE)
  cat("Summary of significant results saved to:", summary_out, "\n")
}

cat("\n========== Analysis Complete ==========\n")
