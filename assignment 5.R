# ============================================================
# Differential Expression Analysis - GSE2034 (Relapse vs Non-Relapse)
# ============================================================

# ---- Load Libraries ----
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hgu133a.db)
library(ggplot2)
library(pheatmap)
library(dplyr)

# ---- Create Results Folder ----
if (!dir.exists("Results")) dir.create("Results")

# ---- 1. Load Expression Data ----
gse <- getGEO("GSE2034", GSEMatrix = TRUE)
gse <- gse[[1]]

exprs_data <- exprs(gse)
targets <- pData(gse)

# Extract group information (Relapse status)
targets$Condition <- ifelse(targets$`disease outcome:ch1` == "Relapse", "Relapse", "NonRelapse")
targets$Condition <- factor(targets$Condition, levels = c("NonRelapse", "Relapse"))

# ---- 2. Probe to Gene Mapping ----
gene_symbols <- mapIds(hgu133a.db,
                       keys = rownames(exprs_data),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Add gene symbols
exprs_data <- as.data.frame(exprs_data)
exprs_data$GeneSymbol <- gene_symbols

# Count duplicates
dup_counts <- table(gene_symbols)
multiple_probes <- sum(dup_counts > 1)
cat("Number of genes with multiple probes:", multiple_probes, "\n")

# Remove probes without gene symbols
exprs_data <- exprs_data[!is.na(exprs_data$GeneSymbol), ]

# Average duplicate probes
exprs_data <- exprs_data %>%
  group_by(GeneSymbol) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  as.data.frame()

# Convert to matrix for limma
rownames(exprs_data) <- exprs_data$GeneSymbol
exprs_matrix <- as.matrix(exprs_data[,-1])

# ---- 3. Differential Gene Expression (Limma) ----
design <- model.matrix(~ 0 + targets$Condition)
colnames(design) <- levels(targets$Condition)
fit <- lmFit(exprs_matrix, design)
contrast.matrix <- makeContrasts(Relapse_vs_NonRelapse = Relapse - NonRelapse, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, coef = "Relapse_vs_NonRelapse", number = Inf, adjust = "fdr")

# ---- 4. Save DEG Results ----
write.csv(deg_results, "Results/DEG_results_complete.csv")

upregulated <- subset(deg_results, logFC > 1 & adj.P.Val < 0.05)
downregulated <- subset(deg_results, logFC < -1 & adj.P.Val < 0.05)

write.csv(upregulated, "Results/DEG_upregulated.csv")
write.csv(downregulated, "Results/DEG_downregulated.csv")

cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# ---- 5. Volcano Plot ----
deg_results$threshold <- as.factor(ifelse(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) >= 1,
                                          ifelse(deg_results$logFC > 1, "Up", "Down"), "Not Sig"))

volcano <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Relapse vs Non-Relapse (GSE2034)",
       x = "log2 Fold Change",
       y = "-log10 Adjusted p-value")

ggsave("Results/VolcanoPlot.png", plot = volcano, width = 7, height = 5)

# ---- 6. Heatmap of Top 25 DEGs ----
top25 <- head(deg_results, 25)
heat_data <- exprs_matrix[rownames(exprs_matrix) %in% rownames(top25), ]

png("Results/Top25_Heatmap.png", width = 800, height = 600)
pheatmap(heat_data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         annotation_col = data.frame(Group = targets$Condition),
         main = "Top 25 Differentially Expressed Genes (GSE2034)")
dev.off()

cat("✅ Analysis complete. Results saved in 'Results/' folder.\n")
