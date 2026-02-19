############################################################
# vc8w RNA-seq analysis: Liver (lvr) (and Small Intestine (si))
############################################################

#### Setup ####

## Core packages
library(GEOquery)
library(DESeq2)
library(gtools)
library(readr)
library(biobroom)
library(dplyr)
library(ggplot2)
library(ggforce)
library(pheatmap)
library(vsn)

## Enrichment & annotation
library(clusterProfiler)
library(AnnotationHub)
library(annotables)
library(org.Mm.eg.db)
library(msigdbr)

## ssGSEA
library(ssGSEA2)
library(tidyverse)

set.seed(8888)

############################################################
#### 1. Load and preprocess featureCounts data ####
############################################################

temp_vc8w <- read.delim(
  "../Data/fc_8w_vc.txt",
  comment.char = "#"
)
colnames(temp_vc8w) <- gsub("\\.LH.*", "", colnames(temp_vc8w))

# Select gene ID column + sample columns (7:30)
vc8w_comb <- temp_vc8w[, c(1, 7:30)]

# Set gene IDs as rownames
vc8w_comb1           <- vc8w_comb
rownames(vc8w_comb1) <- vc8w_comb1[, 1]
vc8w_comb1           <- vc8w_comb1[, 2:ncol(vc8w_comb1)]

# Build metadata
vc8w_comb1_meta <- data.frame(
  colnames(vc8w_comb1),
  do.call(rbind, strsplit(colnames(vc8w_comb1), "\\.")),
  stringsAsFactors = FALSE
)
colnames(vc8w_comb1_meta) <- c("SampleID", "organ", "genotype", "diet", "id")
vc8w_comb1_meta$genotype  <- factor(vc8w_comb1_meta$genotype, levels = c("WT", "vc"))
vc8w_comb1_meta$diet      <- factor(vc8w_comb1_meta$diet,     levels = c("Ctrl", "LP"))
vc8w_comb1_meta$group     <- interaction(vc8w_comb1_meta$diet, vc8w_comb1_meta$genotype)

# Split by organ
lvr_vc8w      <- vc8w_comb1[, grep("lvr", colnames(vc8w_comb1))]
lvr_vc8w_meta <- vc8w_comb1_meta[grep("lvr", colnames(vc8w_comb1)), ]
rownames(lvr_vc8w_meta) <- lvr_vc8w_meta$SampleID

si_vc8w      <- vc8w_comb1[, grep("si", colnames(vc8w_comb1))]
si_vc8w_meta <- vc8w_comb1_meta[grep("si", colnames(vc8w_comb1)), ]
rownames(si_vc8w_meta)  <- si_vc8w_meta$SampleID

############################################################
#### 2. LIVER ####
############################################################

model.matrix(~ group, lvr_vc8w_meta)

dds_lvr_vc8w <- DESeqDataSetFromMatrix(
  countData = as.matrix(lvr_vc8w),
  colData   = as.data.frame(lvr_vc8w_meta),
  design    = ~ group
)
dds_lvr_vc8w <- DESeq(dds_lvr_vc8w)

res.lvr_vc8w          <- results(dds_lvr_vc8w)
t.res.lvr_vc8w        <- tidy.DESeqResults(res.lvr_vc8w)
t.res.lvr_vc8w        <- arrange(t.res.lvr_vc8w, p.adjusted)
resultsNames(dds_lvr_vc8w)

#### 2a. Liver VST & PCA ####

vds_lvr_vc8w <- vst(dds_lvr_vc8w, blind = FALSE)
vsn::meanSdPlot(assay(vds_lvr_vc8w))

plotPCA(vds_lvr_vc8w, intgroup = "group") +
  ggforce::geom_mark_ellipse(aes(color = group))

pca_lvr_vc8w <- plotPCA(vds_lvr_vc8w, intgroup = "group", returnData = TRUE)
percentVar    <- round(100 * attr(pca_lvr_vc8w, "percentVar"))

ggplot(pca_lvr_vc8w, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text(hjust = 0, vjust = 0) +
  ggforce::geom_mark_ellipse(aes(fill = group, color = group))

#### 2b. Liver DEGs (all contrasts) ####

get_DEGs <- function(dds, contrast, alpha = 0.99, lfc = log2(1.4), pval_cutoff = 0.05) {
  res  <- results(dds, alpha = alpha, lfcThreshold = lfc, contrast = contrast)
  degs <- res@rownames[res$pvalue < pval_cutoff & !is.na(res$pvalue)]
  message(contrast[2], " vs ", contrast[3], ": ", length(degs), " DEGs")
  summary(res)
  return(degs)
}

DEG_lvr_vc8w_lpwtvctrlwt  <- get_DEGs(dds_lvr_vc8w, c("group", "LP.WT",   "Ctrl.WT"))
DEG_lvr_vc8w_ctrlvcvctrlwt <- get_DEGs(dds_lvr_vc8w, c("group", "Ctrl.vc", "Ctrl.WT"))
DEG_lvr_vc8w_lpvcvctrlwt  <- get_DEGs(dds_lvr_vc8w, c("group", "LP.vc",   "Ctrl.WT"))
DEG_lvr_vc8w_lpvcvlpwt    <- get_DEGs(dds_lvr_vc8w, c("group", "LP.vc",   "LP.WT"))

saveRDS(DEG_lvr_vc8w_lpwtvctrlwt,   "DEG_lvr_vc8w_lpwtvctrlwt.rds")
saveRDS(DEG_lvr_vc8w_ctrlvcvctrlwt, "DEG_lvr_vc8w_ctrlvcvctrlwt.rds")
saveRDS(DEG_lvr_vc8w_lpvcvctrlwt,   "DEG_lvr_vc8w_lpvcvctrlwt.rds")
saveRDS(DEG_lvr_vc8w_lpvcvlpwt,     "DEG_lvr_vc8w_lpvcvlpwt.rds")

DEG_lvr_vc8w <- unique(c(
  DEG_lvr_vc8w_ctrlvcvctrlwt,
  DEG_lvr_vc8w_lpvcvctrlwt,
  DEG_lvr_vc8w_lpvcvlpwt,
  DEG_lvr_vc8w_lpwtvctrlwt
))
saveRDS(DEG_lvr_vc8w, "DEG_lvr_vc8w.rds")

#### 2c. Liver Heatmap ####

df_col_lvr <- data.frame(cluster = colData(dds_lvr_vc8w)[, "group"])
rownames(df_col_lvr) <- colData(dds_lvr_vc8w)[, "SampleID"]

annot_color_lvr <- list(cluster = c(
  Ctrl.WT = "black",
  LP.WT   = "red",
  Ctrl.vc = "grey",
  LP.vc   = "magenta"
))

# Save to PDF
p_R <- pheatmap(
  assay(vds_lvr_vc8w)[DEG_lvr_vc8w, ],
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col_lvr,
  annotation_colors = annot_color_lvr,
  main              = "DEG_lvr_vc8w",
  cutree_rows       = 3,
  breaks            = seq(-2, 2, length.out = 101),
  color             = colorRampPalette(c("navy", "white", "red"))(100),
  filename          = "HM_lvr_vc8w.pdf"
)
dev.off()

# (no filename)
p_R <- pheatmap(
  assay(vds_lvr_vc8w)[DEG_lvr_vc8w, ],
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col_lvr,
  annotation_colors = annot_color_lvr,
  main              = "DEG_lvr_vc8w",
  cutree_rows       = 3,
  breaks            = seq(-2, 2, length.out = 101),
  color             = colorRampPalette(c("navy", "white", "red"))(100)
)
saveRDS(vds_lvr_vc8w, "lvr_vc8w_vst.rds")

#### 2d. Liver Heatmap Clustering & GO enrichment ####

cl_lvr_vc8w <- cutree(p_R$tree_row, 3)

# Visualize each cluster
for (k in 1:3) {
  pheatmap(
    assay(vds_lvr_vc8w)[names(cl_lvr_vc8w)[cl_lvr_vc8w == k], ],
    scale             = "row",
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = TRUE,
    annotation_col    = df_col_lvr,
    annotation_colors = annot_color_lvr,
    main              = paste0("lvr_vc8w Cluster ", k)
  )
}

# Gene annotation per cluster
lvr_vc8w_cl1 <- inner_join(data.frame(gene = names(cl_lvr_vc8w)[cl_lvr_vc8w == 1]), grcm38, by = c("gene" = "ensgene"))
lvr_vc8w_cl2 <- inner_join(data.frame(gene = names(cl_lvr_vc8w)[cl_lvr_vc8w == 2]), grcm38, by = c("gene" = "ensgene"))
lvr_vc8w_cl3 <- inner_join(data.frame(gene = names(cl_lvr_vc8w)[cl_lvr_vc8w == 3]), grcm38, by = c("gene" = "ensgene"))

# GO enrichment per cluster
go_enrichment_params <- list(OrgDb = org.Mm.eg.db, ont = "BP",
                              pAdjustMethod = "BH", keyType = "ENSEMBL",
                              pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

go_lvr_vc8w_1 <- do.call(enrichGO, c(list(gene = names(cl_lvr_vc8w)[cl_lvr_vc8w == 1]), go_enrichment_params))
dotplot(go_lvr_vc8w_1, showCategory = 20, title = "lvr_vc8w Cluster 1")
write.csv(go_lvr_vc8w_1@result, "lvr_vc8w_GO_1.csv")

go_lvr_vc8w_2 <- do.call(enrichGO, c(list(gene = names(cl_lvr_vc8w)[cl_lvr_vc8w == 2]), go_enrichment_params))
dotplot(go_lvr_vc8w_2, showCategory = 50, title = "lvr_vc8w Cluster 2")
write.csv(go_lvr_vc8w_2@result, "lvr_vc8w_GO_2.csv")

go_lvr_vc8w_3 <- do.call(enrichGO, c(list(gene = names(cl_lvr_vc8w)[cl_lvr_vc8w == 3]), go_enrichment_params))
dotplot(go_lvr_vc8w_3, showCategory = 20, title = "lvr_vc8w Cluster 3")
write.csv(go_lvr_vc8w_3@result, "lvr_vc8w_GO_3.csv")


############################################################
#### 3. SMALL INTESTINE: DESeq2 ####
############################################################

model.matrix(~ group, si_vc8w_meta)

dds_si_vc8w <- DESeqDataSetFromMatrix(
  countData = as.matrix(si_vc8w),
  colData   = as.data.frame(si_vc8w_meta),
  design    = ~ group
)
dds_si_vc8w <- DESeq(dds_si_vc8w)

res.si_vc8w    <- results(dds_si_vc8w)
t.res.si_vc8w  <- tidy.DESeqResults(res.si_vc8w)
t.res.si_vc8w  <- arrange(t.res.si_vc8w, p.adjusted)
resultsNames(dds_si_vc8w)

#### 3a. SI VST & PCA ####

vds_si_vc8w <- vst(dds_si_vc8w, blind = FALSE)
vsn::meanSdPlot(assay(vds_si_vc8w))

plotPCA(vds_si_vc8w, intgroup = "group") +
  ggforce::geom_mark_ellipse(aes(color = group))

pca_si_vc8w <- plotPCA(vds_si_vc8w, intgroup = "group", returnData = TRUE)
percentVar   <- round(100 * attr(pca_si_vc8w, "percentVar"))

ggplot(pca_si_vc8w, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text(hjust = 0, vjust = 0) +
  ggforce::geom_mark_ellipse(aes(fill = group, color = group))

#### 3b. SI DEGs (all contrasts) ####

DEG_si_vc8w_lpwtvctrlwt   <- get_DEGs(dds_si_vc8w, c("group", "LP.WT",   "Ctrl.WT"))
DEG_si_vc8w_ctrlvcvctrlwt <- get_DEGs(dds_si_vc8w, c("group", "Ctrl.vc", "Ctrl.WT"))
DEG_si_vc8w_lpvcvctrlwt   <- get_DEGs(dds_si_vc8w, c("group", "LP.vc",   "Ctrl.WT"))
DEG_si_vc8w_lpvcvlpwt     <- get_DEGs(dds_si_vc8w, c("group", "LP.vc",   "LP.WT"))

saveRDS(DEG_si_vc8w_lpwtvctrlwt,   "DEG_si_vc8w_lpwtvctrlwt.rds")
saveRDS(DEG_si_vc8w_ctrlvcvctrlwt, "DEG_si_vc8w_ctrlvcvctrlwt.rds")
saveRDS(DEG_si_vc8w_lpvcvctrlwt,   "DEG_si_vc8w_lpvcvctrlwt.rds")
saveRDS(DEG_si_vc8w_lpvcvlpwt,     "DEG_si_vc8w_lpvcvlpwt.rds")

DEG_si_vc8w <- unique(c(
  DEG_si_vc8w_ctrlvcvctrlwt,
  DEG_si_vc8w_lpvcvctrlwt,
  DEG_si_vc8w_lpvcvlpwt,
  DEG_si_vc8w_lpwtvctrlwt
))
saveRDS(DEG_si_vc8w, "DEG_si_vc8w.rds")

#### 3c. SI Heatmap ####

df_col_si <- data.frame(cluster = colData(dds_si_vc8w)[, "group"])
rownames(df_col_si) <- colData(dds_si_vc8w)[, "SampleID"]

annot_color_si <- list(cluster = c(
  Ctrl.WT = "black",
  LP.WT   = "red",
  Ctrl.vc = "grey",
  LP.vc   = "magenta"
))

# Save to PDF
p_si <- pheatmap(
  assay(vds_si_vc8w)[DEG_si_vc8w, ],
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col_si,
  annotation_colors = annot_color_si,
  main              = "DEG_si_vc8w",
  cutree_rows       = 4,
  breaks            = seq(-2, 2, length.out = 101),
  color             = colorRampPalette(c("navy", "white", "red"))(100),
  filename          = "HM_si_vc8w.pdf"
)
dev.off()

# Interactive plot (no filename)
p_si <- pheatmap(
  assay(vds_si_vc8w)[DEG_si_vc8w, ],
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col_si,
  annotation_colors = annot_color_si,
  main              = "DEG_si_vc8w",
  cutree_rows       = 4,
  breaks            = seq(-2, 2, length.out = 101),
  color             = colorRampPalette(c("navy", "white", "red"))(100)
)

#### 3d. SI Heatmap Clustering & GO enrichment ####

cl_si_vc8w <- cutree(p_si$tree_row, 4)

# Visualize each cluster
for (k in 1:4) {
  pheatmap(
    assay(vds_si_vc8w)[names(cl_si_vc8w)[cl_si_vc8w == k], ],
    scale             = "row",
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = TRUE,
    annotation_col    = df_col_si,
    annotation_colors = annot_color_si,
    main              = paste0("si_vc8w Cluster ", k)
  )
}

# Gene annotation per cluster
si_vc8w_cl1 <- inner_join(data.frame(gene = names(cl_si_vc8w)[cl_si_vc8w == 1]), grcm38, by = c("gene" = "ensgene"))
si_vc8w_cl2 <- inner_join(data.frame(gene = names(cl_si_vc8w)[cl_si_vc8w == 2]), grcm38, by = c("gene" = "ensgene"))
si_vc8w_cl3 <- inner_join(data.frame(gene = names(cl_si_vc8w)[cl_si_vc8w == 3]), grcm38, by = c("gene" = "ensgene"))
si_vc8w_cl4 <- inner_join(data.frame(gene = names(cl_si_vc8w)[cl_si_vc8w == 4]), grcm38, by = c("gene" = "ensgene"))

# GO enrichment per cluster
go_si_vc8w_1 <- do.call(enrichGO, c(list(gene = names(cl_si_vc8w)[cl_si_vc8w == 1]), go_enrichment_params))
dotplot(go_si_vc8w_1, showCategory = 20, title = "si_vc8w Cluster 1")
write.csv(go_si_vc8w_1@result, "si_vc8w_GO_1.csv")

go_si_vc8w_2 <- do.call(enrichGO, c(list(gene = names(cl_si_vc8w)[cl_si_vc8w == 2]), go_enrichment_params))
dotplot(go_si_vc8w_2, showCategory = 50, title = "si_vc8w Cluster 2")
write.csv(go_si_vc8w_2@result, "si_vc8w_GO_2.csv")

go_si_vc8w_3 <- do.call(enrichGO, c(list(gene = names(cl_si_vc8w)[cl_si_vc8w == 3]), go_enrichment_params))
dotplot(go_si_vc8w_3, showCategory = 20, title = "si_vc8w Cluster 3")
write.csv(go_si_vc8w_3@result, "si_vc8w_GO_3.csv")

go_si_vc8w_4 <- do.call(enrichGO, c(list(gene = names(cl_si_vc8w)[cl_si_vc8w == 4]), go_enrichment_params))
dotplot(go_si_vc8w_4, showCategory = 20, title = "si_vc8w Cluster 4")
write.csv(go_si_vc8w_4@result, "si_vc8w_GO_4.csv")

############################################################
# End of script
############################################################