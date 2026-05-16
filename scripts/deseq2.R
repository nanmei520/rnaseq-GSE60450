# ================================================
# deseq2.R
# 功能：差异表达分析 + 可视化 + 富集分析
# 工具：DESeq2, ggplot2, pheatmap, clusterProfiler
# 数据：GSE60450 | 小鼠乳腺 | chr19 子集
# 比较：Luminal vs Basal
# ================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

setwd("~/rnaseq_practice")

# ════════════════════════════════════════════════
# 第一部分：DESeq2 差异表达分析
# ════════════════════════════════════════════════

# ─── 1. 读取计数矩阵 ────────────────────────────
counts <- read.table("counts/counts_matrix.txt",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t")

# 清理列名：去掉路径和后缀，只保留 SRR ID
colnames(counts) <- gsub("bam/|_Aligned.*", "", colnames(counts))

cat("计数矩阵维度（基因数 × 样本数）：\n")
print(dim(counts))
head(counts)

# ─── 2. 构建样本信息表 ──────────────────────────
sample_info <- data.frame(
    row.names  = c("SRR1552444", "SRR1552445",
                   "SRR1552446", "SRR1552447"),
    condition  = factor(c("Basal", "Basal",
                          "Luminal", "Luminal")),
    sample_name = c("Basal_1", "Basal_2",
                    "Luminal_1", "Luminal_2")
)
print(sample_info)

# ─── 3. 创建 DESeqDataSet 对象 ───────────────────
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sample_info,
    design    = ~ condition
)

# ─── 4. 预过滤：去除低表达基因 ──────────────────
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("过滤后基因数：", nrow(dds), "\n")

# ─── 5. 设置对照组（Basal 为对照）──────────────
dds$condition <- relevel(dds$condition, ref = "Basal")

# ─── 6. 运行 DESeq2 ─────────────────────────────
dds <- DESeq(dds)
cat("Size Factors（测序深度校正系数）：\n")
print(sizeFactors(dds))

# ─── 7. 提取结果（Luminal vs Basal）────────────
res <- results(dds,
               contrast = c("condition", "Luminal", "Basal"),
               alpha    = 0.05)
summary(res)

# ─── 8. LFC 收缩（apeglm 方法）──────────────────
library(apeglm)
res_shrunk <- lfcShrink(dds,
                        coef = "condition_Luminal_vs_Basal",
                        type = "apeglm")

# ─── 9. 保存结果表格 ────────────────────────────
res_df      <- as.data.frame(res_shrunk)
res_df$gene <- rownames(res_df)
res_df      <- res_df[order(res_df$padj), ]

write.csv(res_df,
          file      = "results/DESeq2_results.csv",
          row.names = FALSE)

# 提取显著差异基因（FDR < 0.05，|LFC| > 1）
sig_genes <- res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

cat("显著差异基因数：", nrow(sig_genes), "\n")
head(sig_genes)

# ════════════════════════════════════════════════
# 第二部分：可视化
# ════════════════════════════════════════════════

# ─── PCA 图 ─────────────────────────────────────
rld <- rlog(dds, blind = TRUE)

pca_plot <- plotPCA(rld, intgroup = "condition") +
    theme_bw(base_size = 14) +
    ggtitle("PCA: Sample Clustering (GSE60450 chr19)") +
    scale_color_manual(values = c("Basal"   = "#E74C3C",
                                  "Luminal" = "#3498DB"))

ggsave("results/figures/PCA_plot.pdf",  pca_plot, width = 6, height = 5)
ggsave("results/figures/PCA_plot.png",  pca_plot, width = 6, height = 5, dpi = 300)
print(pca_plot)

# ─── 火山图 ──────────────────────────────────────
volcano <- EnhancedVolcano(res_df,
    lab      = res_df$gene,
    x        = "log2FoldChange",
    y        = "padj",
    xlab     = bquote(~Log[2]~ "Fold Change"),
    title    = "Luminal vs Basal (GSE60450 chr19)",
    pCutoff  = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize   = 4.0,
    col       = c("grey60", "steelblue", "orange", "red3"),
    colAlpha  = 0.7)

ggsave("results/figures/volcano_plot.pdf", volcano, width = 10, height = 8)
ggsave("results/figures/volcano_plot.png", volcano, width = 10, height = 8, dpi = 300)

# ─── 热图 ────────────────────────────────────────
if (nrow(sig_genes) >= 10) {
    top_n   <- min(30, nrow(sig_genes))
    top_genes <- head(sig_genes$gene, top_n)

    mat        <- assay(rld)[top_genes, ]
    mat_scaled <- t(scale(t(mat)))

    annotation_col <- data.frame(
        Condition = sample_info$condition,
        row.names = rownames(sample_info)
    )

    pheatmap(mat_scaled,
             annotation_col = annotation_col,
             color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             show_rownames  = TRUE,
             show_colnames  = TRUE,
             cluster_rows   = TRUE,
             cluster_cols   = TRUE,
             fontsize       = 10,
             main           = paste0("Top ", top_n, " DE Genes: Luminal vs Basal"),
             filename       = "results/figures/heatmap_top30.pdf")
} else {
    cat("显著差异基因数不足 10 个，跳过热图绘制\n")
}

# ─── MA 图 ───────────────────────────────────────
pdf("results/figures/MA_plot.pdf", width = 7, height = 5)
plotMA(res_shrunk, ylim = c(-5, 5),
       main = "MA Plot: LFC Shrinkage (Luminal vs Basal)")
abline(h = c(-1, 1), col = "orange", lty = 2)
dev.off()

# ════════════════════════════════════════════════
# 第三部分：GO / KEGG 富集分析
# ════════════════════════════════════════════════

if (nrow(sig_genes) > 0) {

    # 基因 Symbol → Entrez ID 转换
    entrez_ids <- bitr(
        sig_genes$gene,
        fromType = "ENSEMBL",
        toType   = "ENTREZID",
        OrgDb    = org.Mm.eg.db
    )
    cat("成功转换的基因数：", nrow(entrez_ids), "/", nrow(sig_genes), "\n")

    # ─── GO 富集分析（生物学过程 BP）────────────
    go_result <- enrichGO(
        gene          = entrez_ids$ENTREZID,
        OrgDb         = org.Mm.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )

    if (!is.null(go_result) && nrow(as.data.frame(go_result)) > 0) {
        p_go <- dotplot(go_result, showCategory = 15) +
            ggtitle("GO Biological Process Enrichment")
        ggsave("results/figures/GO_BP_dotplot.pdf", p_go, width = 9, height = 8)
        ggsave("results/figures/GO_BP_dotplot.png", p_go, width = 9, height = 8, dpi = 300)
        write.csv(as.data.frame(go_result), "results/GO_BP_results.csv", row.names = FALSE)
    } else {
        cat("GO 富集分析无显著结果（基因数不足）\n")
    }

    # ─── KEGG 通路富集分析 ───────────────────────
    kegg_result <- enrichKEGG(
        gene          = entrez_ids$ENTREZID,
        organism      = "mmu",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05
    )

    if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
        p_kegg <- barplot(kegg_result, showCategory = 15) +
            ggtitle("KEGG Pathway Enrichment")
        ggsave("results/figures/KEGG_barplot.pdf", p_kegg, width = 9, height = 7)
        ggsave("results/figures/KEGG_barplot.png", p_kegg, width = 9, height = 7, dpi = 300)
        write.csv(as.data.frame(kegg_result), "results/KEGG_results.csv", row.names = FALSE)
    } else {
        cat("KEGG 富集分析无显著结果（基因数不足）\n")
    }

} else {
    cat("无显著差异基因，跳过富集分析\n")
}

cat("\n=== 全部分析完成！结果保存在 results/ 目录 ===\n")
