# ============================================================
# methylGSA Enrichment Analysis - GO / KEGG / Reactome
#
# Input : A CSV file of CpG-level EWAS results containing columns:
#           - unit      : CpG probe ID (e.g. cg00000029)
#           - P.Value   : nominal p-value from the EWAS model
#           - adj.P.Val : BH-adjusted p-value
# Output: GO, KEGG, and Reactome enrichment result CSVs,
#         plus a combined summary of significant pathways.
# ============================================================

# ---- 1. Install required packages ---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Upgrade BiocManager to the version matching your R installation
BiocManager::install(version = "3.23")

# Install methylGSA and its dependencies
BiocManager::install("methylGSA",                                        ask = FALSE)
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",    ask = FALSE)
BiocManager::install("org.Hs.eg.db",                                     ask = FALSE)
BiocManager::install("reactome.db",                                      ask = FALSE)

# ---- 2. Load packages ----------------------------------------
library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ---- 3. Read data --------------------------------------------
# Set input file path and output directory before running
input_file <- "/path/to/your/sig_cpgs_results.csv"
output_dir <- "/path/to/your/output_directory/"

dat <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Inspect structure to confirm column names
head(dat)
cat("Column names:", colnames(dat), "\n")
cat("Total rows:",   nrow(dat),     "\n")

# ---- 4. Extract CpG information ------------------------------
# Expected columns: "unit" (CpG probe ID) and "adj.P.Val" (BH-adjusted p-value)

# All tested CpGs — used as the background set
all_cpg <- dat$unit

# Significant CpGs at FDR < 0.05
sig_cpg <- dat$unit[dat$adj.P.Val < 0.05]

cat("Total CpGs:",                     length(all_cpg), "\n")
cat("Significant CpGs (adj.P < 0.05):", length(sig_cpg), "\n")

# Fall back to nominal p-value threshold if too few significant CpGs are found
if (length(sig_cpg) < 10) {
  cat("Too few significant CpGs; switching to nominal P.Value < 0.05\n")
  sig_cpg <- dat$unit[dat$P.Value < 0.05]
  cat("Significant CpGs after re-filtering:", length(sig_cpg), "\n")
}

# ---- 5. Build named p-value vector (methylGSA input format) --
# methylglm / methylRRA require a named numeric vector of p-values for all CpGs
cpg_pval <- setNames(dat$adj.P.Val, dat$unit)

# ---- 6. GO enrichment analysis -------------------------------
cat("\n========== GO Enrichment ==========\n")

# Primary method: methylglm (logistic regression; recommended for 450K/EPIC data)
# Fallback: methylgometh (hypergeometric test) if methylglm fails
res_GO <- tryCatch({
  methylglm(
    cpg.pval   = cpg_pval,
    array.type = "450K",   # Change to "EPIC" if using the EPIC array
    group      = "all",    # Use all CpGs as background
    GS.type    = "GO",
    minsize    = 5,
    maxsize    = 500
  )
}, error = function(e) {
  cat("methylglm (GO) failed; falling back to methylgometh\n",
      conditionMessage(e), "\n")
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

go_out <- file.path(output_dir, "GO_enrichment_results.csv")
write.csv(res_GO[order(res_GO$pvalue), ], go_out, row.names = FALSE)
cat("GO results saved to:", go_out, "\n")

# ---- 7. KEGG enrichment analysis -----------------------------
cat("\n========== KEGG Enrichment ==========\n")

res_KEGG <- tryCatch({
  methylglm(
    cpg.pval   = cpg_pval,
    array.type = "450K",
    group      = "all",
    GS.type    = "KEGG",
    minsize    = 5,
    maxsize    = 500
  )
}, error = function(e) {
  cat("methylglm (KEGG) failed; falling back to methylgometh\n",
      conditionMessage(e), "\n")
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

# ---- 8. Reactome enrichment analysis -------------------------
cat("\n========== Reactome Enrichment ==========\n")

res_Reactome <- tryCatch({
  methylglm(
    cpg.pval   = cpg_pval,
    array.type = "450K",
    group      = "all",
    GS.type    = "Reactome",
    minsize    = 5,
    maxsize    = 500
  )
}, error = function(e) {
  cat("methylglm (Reactome) failed; falling back to methylgometh\n",
      conditionMessage(e), "\n")
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

# ---- 9. Summarise significant enrichment results (padj < 0.05) ----
cat("\n========== Significant Results Summary (padj < 0.05) ==========\n")

# Helper: filter to significant pathways, handling different column name conventions
# methylGSA may return the adjusted p-value column as "padj", "BH", or "p.adjust"
summarize_sig <- function(res, name) {
  if (is.null(res) || nrow(res) == 0) return(NULL)

  padj_col <- intersect(c("padj", "BH", "p.adjust"), colnames(res))
  if (length(padj_col) == 0) {
    cat(name, ": no adjusted p-value column found; showing entries with pvalue < 0.05\n")
    sig <- res[res$pvalue < 0.05, ]
  } else {
    sig <- res[res[[padj_col[1]]] < 0.05, ]
  }

  cat(name, "significant entries:", nrow(sig), "\n")
  sig$database <- name
  return(sig)
}

sig_GO       <- summarize_sig(res_GO,       "GO")
sig_KEGG     <- summarize_sig(res_KEGG,     "KEGG")
sig_Reactome <- summarize_sig(res_Reactome, "Reactome")

# Combine and save all significant entries across databases
all_sig <- do.call(rbind, Filter(Negate(is.null), list(sig_GO, sig_KEGG, sig_Reactome)))
if (!is.null(all_sig) && nrow(all_sig) > 0) {
  summary_out <- file.path(output_dir, "All_significant_enrichment_summary.csv")
  write.csv(all_sig, summary_out, row.names = FALSE)
  cat("Combined significant results saved to:", summary_out, "\n")
}

cat("\n========== Analysis complete ==========\n")