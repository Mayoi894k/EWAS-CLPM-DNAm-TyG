# ============================================================
# methylGSA 富集分析 - GO / KEGG / Reactome
# 输入: sig_cpgs_p1_adj_age_sex.csv
# ============================================================

# ---- 1. 安装必要的包 ----------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 先升级 BiocManager 到匹配 R 4.6 的版本
BiocManager::install(version = "3.23")

# 然后再安装 methylGSA 及依赖
BiocManager::install("methylGSA", ask = FALSE)
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE)
BiocManager::install("org.Hs.eg.db", ask = FALSE)
BiocManager::install("reactome.db",  ask = FALSE)

# ---- 2. 加载包 -----------------------------------------------
library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ---- 3. 读取数据 ---------------------------------------------
input_file <- "D:/data/twins data/article/修稿/补充分析_富集/sig_cpgs_p1_adj_age_sex.csv"
output_dir <- "D:/data/twins data/article/修稿/补充分析_富集/"

dat <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# 查看数据结构（确认列名）
head(dat)
cat("列名:", colnames(dat), "\n")
cat("总行数:", nrow(dat), "\n")

# ---- 4. 提取 CpG 信息 ----------------------------------------
# 从图片可知：第一列 "unit" 为 CpG ID，"adj.P.Val" 为校正后P值

# 所有检测的 CpG（用于背景）
all_cpg <- dat$unit

# 显著性 CpG（adj.P.Val < 0.05）
sig_cpg <- dat$unit[dat$adj.P.Val < 0.05]

cat("所有CpG数:", length(all_cpg), "\n")
cat("显著CpG数 (adj.P<0.05):", length(sig_cpg), "\n")

# 如果显著CpG太少，尝试用原始P值筛选（可根据需要调整阈值）
if (length(sig_cpg) < 10) {
  cat("显著CpG数量较少，改用原始 P.Value < 0.05 筛选\n")
  sig_cpg <- dat$unit[dat$P.Value < 0.05]
  cat("重新筛选后显著CpG数:", length(sig_cpg), "\n")
}

# ---- 5. 构建 cpg.pval 向量（methylGSA 推荐输入格式）---------
# methylglm / methylRRA 需要一个命名的P值向量（全部CpG）
cpg_pval <- setNames(dat$adj.P.Val, dat$unit)

# ---- 6. GO 富集分析 ------------------------------------------
cat("\n========== GO 富集分析 ==========\n")

## 方法一：methylglm（逻辑回归，推荐）
res_GO <- tryCatch({
  methylglm(
    cpg.pval    = cpg_pval,
    array.type  = "450K",       # 若为EPIC芯片请改为 "EPIC"
    group       = "all",        # 使用全部CpG作为背景
    GS.type     = "GO",
    minsize     = 5,
    maxsize     = 500
  )
}, error = function(e) {
  cat("methylglm GO 出错，改用 methylgometh\n", conditionMessage(e), "\n")
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

cat("GO 富集结果行数:", nrow(res_GO), "\n")
print(head(res_GO[order(res_GO$pvalue), ], 10))

# 保存 GO 结果
go_out <- file.path(output_dir, "GO_enrichment_results.csv")
write.csv(res_GO[order(res_GO$pvalue), ], go_out, row.names = FALSE)
cat("GO 结果已保存至:", go_out, "\n")

# ---- 7. KEGG 富集分析 ----------------------------------------
cat("\n========== KEGG 富集分析 ==========\n")

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
  cat("methylglm KEGG 出错，改用 methylgometh\n", conditionMessage(e), "\n")
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

cat("KEGG 富集结果行数:", nrow(res_KEGG), "\n")
print(head(res_KEGG[order(res_KEGG$pvalue), ], 10))

kegg_out <- file.path(output_dir, "KEGG_enrichment_results.csv")
write.csv(res_KEGG[order(res_KEGG$pvalue), ], kegg_out, row.names = FALSE)
cat("KEGG 结果已保存至:", kegg_out, "\n")

# ---- 8. Reactome 富集分析 ------------------------------------
cat("\n========== Reactome 富集分析 ==========\n")

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
  cat("methylglm Reactome 出错，改用 methylgometh\n", conditionMessage(e), "\n")
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

cat("Reactome 富集结果行数:", nrow(res_Reactome), "\n")
print(head(res_Reactome[order(res_Reactome$pvalue), ], 10))

reactome_out <- file.path(output_dir, "Reactome_enrichment_results.csv")
write.csv(res_Reactome[order(res_Reactome$pvalue), ], reactome_out, row.names = FALSE)
cat("Reactome 结果已保存至:", reactome_out, "\n")

# ---- 9. 汇总显著富集结果 (padj < 0.05) ----------------------
cat("\n========== 汇总显著结果 (padj < 0.05) ==========\n")

summarize_sig <- function(res, name) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  # methylGSA 结果列名可能是 padj 或 BH
  padj_col <- intersect(c("padj", "BH", "p.adjust"), colnames(res))
  if (length(padj_col) == 0) {
    cat(name, ": 未找到校正P值列，显示P值<0.05的结果\n")
    sig <- res[res$pvalue < 0.05, ]
  } else {
    sig <- res[res[[padj_col[1]]] < 0.05, ]
  }
  cat(name, "显著富集条目数:", nrow(sig), "\n")
  sig$database <- name
  return(sig)
}

sig_GO       <- summarize_sig(res_GO,       "GO")
sig_KEGG     <- summarize_sig(res_KEGG,     "KEGG")
sig_Reactome <- summarize_sig(res_Reactome, "Reactome")

# 合并保存
all_sig <- do.call(rbind, Filter(Negate(is.null), list(sig_GO, sig_KEGG, sig_Reactome)))
if (!is.null(all_sig) && nrow(all_sig) > 0) {
  summary_out <- file.path(output_dir, "All_significant_enrichment_summary.csv")
  write.csv(all_sig, summary_out, row.names = FALSE)
  cat("汇总显著结果已保存至:", summary_out, "\n")
}

cat("\n========== 分析完成 ==========\n")