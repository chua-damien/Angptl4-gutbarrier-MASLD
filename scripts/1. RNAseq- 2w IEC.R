############################################################
# IEC RNA-seq Angptl4 Vil-cre analysis 
############################################################

#### Set up ####

## Core packages for DESeq2-based RNA-seq and plotting
library(GEOquery)
library(DESeq2)
library(gtools)
library(biobroom)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(sva)
library(VennDiagram)

## Enrichment / annotation packages
library(clusterProfiler)
library(AnnotationHub)
library(annotables)
library(org.Mm.eg.db)
library(DOSE)
library(msigdbr)
library(enrichplot)
library(ggnewscale)
library(ggrepel)
library(plyr)
# library(ashr)

set.seed(8888)

############################################################
#### PART 1: Data loading, combining & DESeq2 ####
############################################################

#### Loading datasets ####

# LIDPAD vs Control
temp_LPvCtrl <- read.csv("./Datasets/ZS_LPvCtrl_Fcount/IEC_LPvCtrl_Fcount.csv", row.names = 1)
colnames(temp_LPvCtrl) <- temp_LPvCtrl[1, ]
temp_LPvCtrl <- temp_LPvCtrl[2:nrow(temp_LPvCtrl), ]
colnames(temp_LPvCtrl) <- sub("\\\\LE.*", "", colnames(temp_LPvCtrl))
LPvCtrl <- temp_LPvCtrl[, 6:ncol(temp_LPvCtrl)]

LPvCtrl_cond <- data.frame(
  colnames(LPvCtrl),
  conds = factor(c(rep("Ctrl", 4), rep("LP", 4)), levels = c("Ctrl", "LP"))
)

# Angptl4 Vil-cre IEC Tx
temp_fiafvilcre <- read.table("./Datasets/Angptl4-Vil-cre_IEC_Fcount/Fcount_Vilcre_MP.txt", row.names = 1)
colnames(temp_fiafvilcre) <- temp_fiafvilcre[1, ]
temp_fiafvilcre <- temp_fiafvilcre[2:nrow(temp_fiafvilcre), ]

# Arrange dataset
colnames(temp_fiafvilcre) <- sub("\\\\-.*", "", colnames(temp_fiafvilcre))
temp_fiafvilcre <- temp_fiafvilcre[, c(
  "Chr", "Start", "End", "Strand", "Length",
  "MP14", "MP16a",
  "LG42", "LG48", "MP7", "MP8", # LP
  "LG44", "LG52",               # Vil-cre Ctrl
  "LG40", "LG41", "LG49", "LG51" # Vil-cre LP
)]

# Remove chr, start, end, strand, length
fiafvilcre <- temp_fiafvilcre[, 6:ncol(temp_fiafvilcre)]

fiafvilcre_cond <- factor(c(
  rep("Ctrl", 2), rep("LP", 4),
  rep("Fiafvc-Ctrl", 2), rep("Fiafvc-LP", 4)
))

# Mega combine mode: merge feature-count matrices

temp_fiafvilcre1 <- cbind(rownames(temp_fiafvilcre), temp_fiafvilcre)
colnames(temp_fiafvilcre1)[1] <- "gene.id"

LPvCtrl1 <- cbind(rownames(LPvCtrl), LPvCtrl)
colnames(LPvCtrl1)[1] <- "gene.id"

temp_fiafvilcre_comb <- merge(
  temp_fiafvilcre1,
  LPvCtrl1,
  by = 1,
  all = FALSE
)

# Arrange dataset according to conditions
temp_fiafvilcre_comb <- temp_fiafvilcre_comb[, c(
  "gene.id", "Chr", "Start", "End", "Strand", "Length",
  "MP14", "MP16a",
  "1r2w-c-1-", "1r2w-c-2-", "1r2w-c-3-", "1r2w-c-4-", # Ctrl
  "LG42", "LG48", "MP7", "MP8", "1r2w-lp1-", "1r2w-lp2-", "1r2w-lp3-", "1r2w-lp4-", # LP
  "LG44", "LG52", # Vil-cre Ctrl
  "LG40", "LG49", "LG51" # Vil-cre LP
)]

fiafvc_comb <- temp_fiafvilcre_comb[, c(1, 7:ncol(temp_fiafvilcre_comb))]
rownames(fiafvc_comb) <- fiafvc_comb[, 1]
fiafvc_comb <- fiafvc_comb[, 2:ncol(fiafvc_comb)]

fiafvc_comb_cond <- data.frame(
  colnames(fiafvc_comb),
  genotype = factor(c(
    rep("flfl", 6), rep("flfl", 8),
    rep("Fiafvc", 2), rep("Fiafvc", 3)
  ), levels = c("flfl", "Fiafvc")),
  conds = factor(c(
    rep("Ctrl", 6), rep("LP", 8),
    rep("Ctrl", 2), rep("LP", 3)
  ), levels = c("Ctrl", "LP")),
  batch = factor(c(
    rep(1, 2), rep(2, 4),
    rep(3, 2), rep(1, 2), rep(2, 4),
    rep(3, 2), rep(3, 3)
  )),
  cluster = factor(c(
    rep("flfl-Ctrl", 6), rep("flfl-LP", 8),
    rep("fiafvc-Ctrl", 2), rep("fiafvc-LP", 3)
  ), levels = c("flfl-Ctrl", "flfl-LP", "fiafvc-Ctrl", "fiafvc-LP"))
)

#### Batch correction with ComBat-seq ####

# Convert counts to numeric
fiafvc_comb1 <- dplyr::mutate_all(fiafvc_comb, function(x) as.numeric(as.character(x)))

covar_mat <- data.frame(
  genotype = fiafvc_comb_cond$genotype,
  conds    = fiafvc_comb_cond$conds,
  cluster  = fiafvc_comb_cond$cluster
)
covar_mat <- as.matrix(covar_mat)
rownames(covar_mat) <- fiafvc_comb_cond$colnames.fiafvc_comb.

fiafvc_comb_cbs <- ComBat_seq(
  as.matrix(fiafvc_comb1),
  batch = fiafvc_comb_cond$batch,
  group = covar_mat
)

#### DESeq2 on ComBat-seq corrected counts ####

model.matrix(~ batch + genotype + conds + genotype:conds, fiafvc_comb_cond)

dds_fiafvc_comb_cbs <- DESeqDataSetFromMatrix(
  countData = as.matrix(fiafvc_comb_cbs),
  colData   = as.data.frame(fiafvc_comb_cond),
  design    = ~ cluster
)

dds_fiafvc_comb_cbs <- DESeq(dds_fiafvc_comb_cbs)
resultsNames(dds_fiafvc_comb_cbs)

#### Variance-stabilising transformation & PCA ####

vds_fiafvc_comb_cbs <- vst(dds_fiafvc_comb_cbs, blind = FALSE)
vsn::meanSdPlot(assay(vds_fiafvc_comb_cbs))

plotPCA(vds_fiafvc_comb_cbs, intgroup = c("batch"))
plotPCA(vds_fiafvc_comb_cbs, intgroup = c("genotype"))
plotPCA(vds_fiafvc_comb_cbs, intgroup = c("conds", "genotype", "batch")) +
  ggforce::geom_mark_ellipse(aes(
    fill  = conds,
    color = conds,
    shape = genotype
  ))

pca_fiafvc_comb_cbs <- plotPCA(
  vds_fiafvc_comb_cbs,
  intgroup   = c("conds", "genotype", "batch", "cluster"),
  returnData = TRUE
)

percentVar <- round(100 * attr(pca_fiafvc_comb_cbs, "percentVar"))
ggplot(
  pca_fiafvc_comb_cbs,
  aes(PC1, PC2, color = conds, label = name, shape = genotype)
) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text(hjust = 0, vjust = 0) +
  ggforce::geom_mark_ellipse(aes(
    fill  = cluster,
    color = cluster,
    label = cluster
  ))

#### Heatmap with DEG_LPvCtrl ####

library(pheatmap)
DEG_LPvCtrl <- readRDS("DEG_LPvCtrl.rds")

df_col <- data.frame(cluster = colData(dds_fiafvc_comb_cbs)[, c("cluster")])
rownames(df_col) <- colData(dds_fiafvc_comb_cbs)[, c("colnames.fiafvc_comb.")]

annot_color <- c(
  rep("black", 6),
  rep("red", 8),
  rep("grey", 2),
  rep("magenta", 3)
)
names(annot_color) <- df_col$cluster
annot_color <- list(cluster = annot_color)

p_fiafvc_cb <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[DEG_LPvCtrl, ],
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)

p_fiafvc_cb <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[
    DEG_LPvCtrl,
    which(colData(vds_fiafvc_comb_cbs)$cluster %in% c("flfl-Ctrl", "flfl-LP"))
  ],
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = TRUE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)

#### DEG lists for different contrasts ####

resultsNames(dds_fiafvc_comb_cbs)

## flfl LP vs Ctrl
res.fiafvc_comb_cbs <- results(
  dds_fiafvc_comb_cbs,
  contrast = c("cluster", "flfl-LP", "flfl-Ctrl")
)
t.res.fiafvc_comb_cbs <- tidy.DESeqResults(res.fiafvc_comb_cbs)
t.res.fiafvc_comb_cbs <- arrange(t.res.fiafvc_comb_cbs, p.adjusted)

res.fiafvc_comb_cbs_filt <- results(
  dds_fiafvc_comb_cbs,
  alpha    = 0.05,
  contrast = c("cluster", "flfl-LP", "flfl-Ctrl")
)
DEG_comb_CBS_flfl_LPvCtrl <- res.fiafvc_comb_cbs_filt@rownames[
  res.fiafvc_comb_cbs_filt$padj < 0.05 & !is.na(res.fiafvc_comb_cbs_filt$padj)
]
saveRDS(DEG_comb_CBS_flfl_LPvCtrl, "DEG_comb_CBS_flfl_LPvCtrl.rds")

## fiafvc LP vs Ctrl
res.fiafvc_comb_cbs_filt <- results(
  dds_fiafvc_comb_cbs,
  alpha    = 0.05,
  contrast = c("cluster", "fiafvc-LP", "fiafvc-Ctrl")
)
DEG_comb_CBS_fiafvc_LPvCtrl <- res.fiafvc_comb_cbs_filt@rownames[
  res.fiafvc_comb_cbs_filt$padj < 0.05 & !is.na(res.fiafvc_comb_cbs_filt$padj)
]
saveRDS(DEG_comb_CBS_fiafvc_LPvCtrl, "DEG_comb_CBS_fiafvc_LPvCtrl.rds")

## Ctrl fiafvc vs flfl
res.fiafvc_comb_cbs_filt <- results(
  dds_fiafvc_comb_cbs,
  alpha    = 0.05,
  contrast = c("cluster", "fiafvc-Ctrl", "flfl-Ctrl")
)
DEG_comb_CBS_Ctrl_fiafvcvflfl <- res.fiafvc_comb_cbs_filt@rownames[
  res.fiafvc_comb_cbs_filt$padj < 0.05 & !is.na(res.fiafvc_comb_cbs_filt$padj)
]

## LP fiafvc vs flfl
res.fiafvc_comb_cbs_filt <- results(
  dds_fiafvc_comb_cbs,
  alpha    = 0.05,
  contrast = c("cluster", "fiafvc-LP", "flfl-LP")
)
DEG_comb_CBS_LP_fiafvcvflfl <- res.fiafvc_comb_cbs_filt@rownames[
  res.fiafvc_comb_cbs_filt$padj < 0.05 & !is.na(res.fiafvc_comb_cbs_filt$padj)
]

#### Remove variant genes for genotype contrasts ####

df_col <- data.frame(cluster = colData(dds_fiafvc_comb_cbs)[, c("cluster")])
rownames(df_col) <- colData(dds_fiafvc_comb_cbs)[, c("colnames.fiafvc_comb.")]
annot_color <- c(
  rep("black", 6),
  rep("red", 8),
  rep("grey", 2),
  rep("magenta", 3)
)
names(annot_color) <- df_col$cluster
annot_color <- list(cluster = annot_color)

## Fiafvc Ctrl vs flfl Ctrl
set.seed(8888)
p_comb_cbs <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[DEG_comb_CBS_Ctrl_fiafvcvflfl, ],
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)

cl_comb_cbs_Ctrl_fiafvcvflfl <- cutree(p_comb_cbs_Ctrl_fiafvcvflfl$tree_row, 9)
variant2 <- DEG_comb_CBS_Ctrl_fiafvcvflfl[
  DEG_comb_CBS_Ctrl_fiafvcvflfl %in%
    names(cl_comb_cbs_Ctrl_fiafvcvflfl)[cl_comb_cbs_Ctrl_fiafvcvflfl %in% c(1, 3)]
]
DEG_comb_CBS_Ctrl_fiafvcvflfl <- names(cl_comb_cbs_Ctrl_fiafvcvflfl)[
  cl_comb_cbs_Ctrl_fiafvcvflfl %in% c(2, 4, 5, 6, 7, 8, 9)
]
saveRDS(DEG_comb_CBS_Ctrl_fiafvcvflfl, "DEG_comb_CBS_Ctrl_fiafvcvflfl.rds")

## Fiafvc LP vs flfl LP
set.seed(8888)
p_comb_cbs_LP_fiafvcvflfl <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[DEG_comb_CBS_LP_fiafvcvflfl, ],
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cutree_rows       = 20,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)
cl_comb_cbs_LP_fiafvcvflfl <- cutree(p_comb_cbs_LP_fiafvcvflfl$tree_row, 20)
variant3 <- DEG_comb_CBS_LP_fiafvcvflfl[
  DEG_comb_CBS_LP_fiafvcvflfl %in%
    names(cl_comb_cbs_LP_fiafvcvflfl)[cl_comb_cbs_LP_fiafvcvflfl %in% c(16, 19)]
]
DEG_comb_CBS_LP_fiafvcvflfl <- names(cl_comb_cbs_LP_fiafvcvflfl)[
  !cl_comb_cbs_LP_fiafvcvflfl %in% c(16, 19)
]
saveRDS(DEG_comb_CBS_LP_fiafvcvflfl, "DEG_comb_CBS_LP_fiafvcvflfl.rds")

## Variant DEG list
variant1 <- DEG_comb_CBS_fiafvc_LPvCtrl[
  DEG_comb_CBS_fiafvc_LPvCtrl %in%
    names(cl_comb_cbs_LP_fiafvcvflfl)[cl_comb_cbs_LP_fiafvcvflfl %in% c(1, 3)]
]  # as per original logic (if needed)
DEG_variant <- union(variant1, union(variant2, variant3))
saveRDS(DEG_variant, "DEG_variant.rds")

############################################################
#### PART 2: Venn, diet/genotype/effect genes ####
############################################################

# Continue from previous section: dds_fiafvc_comb_cbs, vds_fiafvc_comb_cbs,
# DEG_comb_CBS_* objects and DEG_variant are already defined.

#### Venn diagrams ####

venn.diagram(
  x = list(
    DEG_comb_CBS_flfl_LPvCtrl,
    DEG_comb_CBS_fiafvc_LPvCtrl,
    DEG_comb_CBS_Ctrl_fiafvcvflfl,
    DEG_comb_CBS_LP_fiafvcvflfl
  ),
  category.names = c("flfl: LPvCtrl", "Fiafvc: LPvCtrl",
                     "Ctrl: FiafvcvFlfc", "LP: FiafvcvFlfc"),
  filename = "all_DEGs_venn_check.png",
  output   = TRUE
)

venn.diagram(
  x = list(DEG_comb_CBS_flfl_LPvCtrl, DEG_comb_CBS_fiafvc_LPvCtrl),
  category.names = c("flfl", "fiafvc"),
  filename       = "Venn-angptl4-dep-diet-DEG.png",
  output         = TRUE,
  imagetype      = "png",
  height         = 1000,
  width          = 1000,
  resolution     = 300,
  cex            = 0.5,
  fontfamily     = "sans",
  cat.cex        = 0.5,
  cat.default.pos = "outer",
  cat.pos         = c(0, 10),
  cat.dist        = c(0.05, 0.05)
)

#### DIET ONLY EFFECT ####

diet.only.genes <- intersect(DEG_comb_CBS_flfl_LPvCtrl, DEG_comb_CBS_fiafvc_LPvCtrl)
diet.only.genes <- data.frame(gene = diet.only.genes)
diet.only.genes <- inner_join(diet.only.genes, grcm38, by = c("gene" = "ensgene"))
write.csv(diet.only.genes, file = "diet_only_7genes.csv")

df_col <- data.frame(cluster = colData(dds_fiafvc_comb_cbs)[, c("cluster")])
rownames(df_col) <- colData(dds_fiafvc_comb_cbs)[, c("colnames.fiafvc_comb.")]
annot_color <- c(
  rep("black", 6),
  rep("red", 8),
  rep("grey", 2),
  rep("magenta", 3)
)
names(annot_color) <- df_col$cluster
annot_color <- list(cluster = annot_color)

set.seed(8888)
hm_plot <- assay(vds_fiafvc_comb_cbs)[
  diet.only.genes$gene,
  c(
    "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
    "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
    "LG44", "LG52",
    "LG49", "LG51"
  )
]
rownames(hm_plot) <- diet.only.genes$symbol

pheatmap(
  hm_plot,
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = TRUE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)

GO <- enrichGO(
  diet.only.genes$gene,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.3,
  qvalueCutoff  = 0.3,
  readable      = TRUE
)
dotplot(GO)


#### Angptl4 differential DIET EFFECT (effect genes) ####

genotype.diet.genes <- union(genotype.only.genes$gene, diet.only.genes$gene)
effect.flfl.genes   <- setdiff(DEG_comb_CBS_flfl_LPvCtrl, genotype.diet.genes)
effect.fiafvc.genes <- setdiff(DEG_comb_CBS_fiafvc_LPvCtrl, genotype.diet.genes)

res.flfl <- results(
  dds_fiafvc_comb_cbs,
  contrast = c("cluster", "flfl-LP", "flfl-Ctrl")
)
t.res.flfl <- tidy.DESeqResults(res.flfl)
t.res.flfl <- arrange(t.res.flfl, p.adjusted)
t.res.flfl <- inner_join(t.res.flfl, grcm38, by = c("gene" = "ensgene"))

flfl_volp <- data.frame(
  ens   = t.res.flfl$gene,
  gene  = t.res.flfl$symbol,
  log2FC = t.res.flfl$estimate,
  p.adj  = t.res.flfl$p.adjusted
)
flfl_volp$deg <- "ns"
flfl_volp$deg[flfl_volp$ens %in% DEG_comb_CBS_flfl_LPvCtrl] <- "sig"
flfl_volp$deg[flfl_volp$ens %in% effect.flfl.genes]        <- "flfl-spec"
flfl_volp$deg[abs(flfl_volp$log2FC) < log2(1.5)]           <- "ns"
flfl_volp <- flfl_volp[!flfl_volp$ens %in% DEG_variant, ]

res.fiafvc <- results(
  dds_fiafvc_comb_cbs,
  contrast = c("cluster", "fiafvc-LP", "fiafvc-Ctrl")
)
t.res.fiafvc <- tidy.DESeqResults(res.fiafvc)
t.res.fiafvc <- arrange(t.res.fiafvc, p.adjusted)
t.res.fiafvc <- inner_join(t.res.fiafvc, grcm38, by = c("gene" = "ensgene"))

fiafvc_volp <- data.frame(
  ens   = t.res.fiafvc$gene,
  gene  = t.res.fiafvc$symbol,
  log2FC = t.res.fiafvc$estimate,
  p.adj  = t.res.fiafvc$p.adjusted
)
fiafvc_volp$deg <- "ns"
fiafvc_volp$deg[fiafvc_volp$ens %in% DEG_comb_CBS_fiafvc_LPvCtrl] <- "sig"
fiafvc_volp$deg[fiafvc_volp$ens %in% effect.fiafvc.genes]         <- "fiafvc-spec"
fiafvc_volp$deg[abs(fiafvc_volp$log2FC) < log2(1.5)]              <- "ns"
fiafvc_volp <- fiafvc_volp[!fiafvc_volp$ens %in% DEG_variant, ]

ggplot(flfl_volp[, 2:5], aes(x = log2FC, y = -log10(p.adj), col = deg)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(log2(2/3), log2(1.5)), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlim(-10, 20)

ggplot(fiafvc_volp[, 2:5], aes(x = log2FC, y = -log10(p.adj), col = deg)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(log2(2/3), log2(1.5)), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlim(-10, 20)

#### Effect genes heatmap and GO by effect clusters ####

effect.genes <- data.frame(gene = union(effect.flfl.genes, effect.fiafvc.genes))
effect.genes <- inner_join(effect.genes, grcm38, by = c("gene" = "ensgene"))
write.csv(effect.genes, file = "effect_genes.csv")

df_col <- data.frame(cluster = colData(dds_fiafvc_comb_cbs)[, c("cluster")])
rownames(df_col) <- colData(dds_fiafvc_comb_cbs)[, c("colnames.fiafvc_comb.")]
annot_color <- c(
  rep("black", 6),
  rep("red", 8),
  rep("grey", 2),
  rep("magenta", 3)
)
names(annot_color) <- df_col$cluster
annot_color <- list(cluster = annot_color)

set.seed(8888)
hm_plot <- assay(vds_fiafvc_comb_cbs)[
  effect.genes$gene,
  c(
    "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
    "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
    "LG44", "LG52",
    "LG49", "LG51"
  )
]
rownames(hm_plot) <- effect.genes$symbol

cl_effect_only <- pheatmap(
  hm_plot,
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)
rownames(hm_plot) <- effect.genes$gene

cl_effect_only <- pheatmap(
  hm_plot,
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  cutree_rows       = 5,
  annotation_col    = df_col,
  annotation_colors = annot_color
)
cl_effect_only <- cutree(cl_effect_only$tree_row, 5)
table(cl_effect_only)

pheatmap(
  hm_plot,
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color
)

## Effect cluster 1
HM_1 <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[
    names(cl_effect_only)[cl_effect_only == 1],
    c(
      "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
      "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
      "LG44", "LG52",
      "LG49", "LG51"
    )
  ],
  scale             = "row",
  cluster_rows      = FALSE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color,
  main              = "effect_cluster_1"
)

GO_1 <- enrichGO(
  effect.genes$gene[cl_effect_only %in% 1],
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.3,
  qvalueCutoff  = 0.3,
  readable      = TRUE
)
write.csv(GO_1@result, "effect_cluster_1.csv")
dotplot(GO_1, showCategory = 15)

## Effect cluster 2
HM_2 <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[
    names(cl_effect_only)[cl_effect_only == 2],
    c(
      "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
      "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
      "LG44", "LG52",
      "LG49", "LG51"
    )
  ],
  scale             = "row",
  cluster_rows      = TRUE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color,
  main              = "effect_cluster_2"
)

GO_2 <- enrichGO(
  effect.genes$gene[cl_effect_only %in% 2],
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.3,
  qvalueCutoff  = 0.3,
  readable      = TRUE
)
write.csv(GO_2@result, "effect_cluster_2.csv")
dotplot(GO_2, showCategory = 15)

## Effect cluster 3
HM_3 <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[
    names(cl_effect_only)[cl_effect_only == 3],
    c(
      "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
      "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
      "LG44", "LG52",
      "LG49", "LG51"
    )
  ],
  scale             = "row",
  cluster_rows      = FALSE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color,
  main              = "effect_cluster_3"
)

GO_3 <- enrichGO(
  effect.genes$gene[cl_effect_only %in% 3],
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.3,
  qvalueCutoff  = 0.3,
  readable      = TRUE
)
write.csv(GO_3@result, "effect_cluster_3.csv")
dotplot(GO_3, showCategory = 15)

## Effect cluster 4
HM_4 <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[
    names(cl_effect_only)[cl_effect_only == 4],
    c(
      "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
      "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
      "LG44", "LG52",
      "LG49", "LG51"
    )
  ],
  scale             = "row",
  cluster_rows      = FALSE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color,
  main              = "effect_cluster_4"
)

GO_4 <- enrichGO(
  effect.genes$gene[cl_effect_only %in% 4],
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  readable      = TRUE
)
write.csv(GO_4@result, "effect_cluster_4.csv")
dotplot(GO_4, showCategory = 15)

## Effect cluster 5
HM_5 <- pheatmap(
  assay(vds_fiafvc_comb_cbs)[
    names(cl_effect_only)[cl_effect_only == 5],
    c(
      "MP16a", "1r2w-c-3-", "1r2w-c-1-", "MP14", "1r2w-c-4-",
      "1r2w-lp1-", "MP7", "1r2w-lp4-", "LG42", "MP8",
      "LG44", "LG52",
      "LG49", "LG51"
    )
  ],
  scale             = "row",
  cluster_rows      = FALSE,
  show_rownames     = FALSE,
  cluster_cols      = FALSE,
  show_colnames     = TRUE,
  annotation_col    = df_col,
  annotation_colors = annot_color,
  main              = "effect_cluster_5"
)

GO_5 <- enrichGO(
  effect.genes$gene[cl_effect_only %in% 5],
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  keyType       = "ENSEMBL",
  pvalueCutoff  = 0.3,
  qvalueCutoff  = 0.3,
  readable      = TRUE
)
write.csv(GO_5@result, "effect_cluster_5.csv")
dotplot(GO_5, showCategory = 15)

############################################################
#### PART 3: Detailed volcano, GSEA & plotting ####
############################################################

# Note: DEG_* objects, DEG_variant, vds_fiafvc_comb_cbs, dds_fiafvc_comb_cbs
# are available from previous sections.

#### Volcano – shared and effect genes with manual highlights ####

diet.only.genes <- intersect(DEG_comb_CBS_flfl_LPvCtrl, DEG_comb_CBS_fiafvc_LPvCtrl)
diet.only.genes <- data.frame(gene = diet.only.genes)
diet.only.genes <- inner_join(diet.only.genes, grcm38, by = c("gene" = "ensgene"))

genotype.only.genes <- union(
  intersect(DEG_comb_CBS_Ctrl_fiafvcvflfl, DEG_comb_CBS_LP_fiafvcvflfl),
  union(
    setdiff(
      DEG_comb_CBS_Ctrl_fiafvcvflfl,
      union(DEG_comb_CBS_fiafvc_LPvCtrl,
            union(DEG_comb_CBS_flfl_LPvCtrl, DEG_comb_CBS_LP_fiafvcvflfl))
    ),
    setdiff(
      DEG_comb_CBS_LP_fiafvcvflfl,
      union(DEG_comb_CBS_fiafvc_LPvCtrl,
            union(DEG_comb_CBS_flfl_LPvCtrl, DEG_comb_CBS_Ctrl_fiafvcvflfl))
    )
  )
)
genotype.only.genes <- data.frame(gene = genotype.only.genes)
genotype.only.genes <- inner_join(genotype.only.genes, grcm38, by = c("gene" = "ensgene"))

genotype.diet.genes <- union(genotype.only.genes$gene, diet.only.genes$gene)
effect.flfl.genes   <- setdiff(DEG_comb_CBS_flfl_LPvCtrl, genotype.diet.genes)
effect.fiafvc.genes <- setdiff(DEG_comb_CBS_fiafvc_LPvCtrl, genotype.diet.genes)

res.flfl <- results(
  dds_fiafvc_comb_cbs,
  contrast = c("cluster", "flfl-LP", "flfl-Ctrl")
)
t.res.flfl <- tidy.DESeqResults(res.flfl)
t.res.flfl <- arrange(t.res.flfl, p.adjusted)
t.res.flfl <- inner_join(t.res.flfl, grcm38, by = c("gene" = "ensgene"))

flfl_volp <- data.frame(
  ens   = t.res.flfl$gene,
  gene  = t.res.flfl$symbol,
  log2FC = t.res.flfl$estimate,
  p.adj  = t.res.flfl$p.adjusted
)
flfl_volp$deg <- "ns"
flfl_volp$deg[flfl_volp$ens %in% DEG_comb_CBS_flfl_LPvCtrl] <- "sig"
flfl_volp$deg[flfl_volp$gene %in% c("Tfrc")]                 <- "up.down"
flfl_volp$deg[flfl_volp$gene %in% c("Glul", "Htra1", "Timp2")] <- "up.up"
flfl_volp$deg[flfl_volp$gene %in% c("Gm15501", "Npdc1", "Cpne8")] <- "down.up"
flfl_volp$deg[abs(flfl_volp$log2FC) < log2(1.5)]             <- "ns"
flfl_volp <- flfl_volp[!flfl_volp$ens %in% DEG_variant, ]

flfl_volp$shared <- NA
flfl_volp$shared[flfl_volp$gene %in% c("Gm15501")] <- "Gm15501"
flfl_volp$shared[flfl_volp$gene %in% c("Htra1")]   <- "Htra1"
flfl_volp$shared[flfl_volp$gene %in% c("Cpne8")]   <- "Cpne8"
flfl_volp$shared[abs(flfl_volp$log2FC) < log2(1.5)] <- NA

flfl_volp_fr  <- which(flfl_volp$deg %in% c("up.up", "down.up"))
flfl_volp_beh <- setdiff(1:nrow(flfl_volp), flfl_volp_fr)
flfl_volp     <- flfl_volp[c(flfl_volp_beh, flfl_volp_fr), ]

res.fiafvc <- results(
  dds_fiafvc_comb_cbs,
  contrast = c("cluster", "fiafvc-LP", "fiafvc-Ctrl")
)
t.res.fiafvc <- tidy.DESeqResults(res.fiafvc)
t.res.fiafvc <- arrange(t.res.fiafvc, p.adjusted)
t.res.fiafvc <- inner_join(t.res.fiafvc, grcm38, by = c("gene" = "ensgene"))

fiafvc_volp <- data.frame(
  ens   = t.res.fiafvc$gene,
  gene  = t.res.fiafvc$symbol,
  log2FC = t.res.fiafvc$estimate,
  p.adj  = t.res.fiafvc$p.adjusted
)
fiafvc_volp$deg <- "ns"
fiafvc_volp$deg[fiafvc_volp$ens %in% DEG_comb_CBS_fiafvc_LPvCtrl] <- "sig"
fiafvc_volp$deg[fiafvc_volp$gene %in% c("Htra1")]                 <- "up.up"
fiafvc_volp$deg[fiafvc_volp$gene %in% c("Gm15501", "Cpne8")]      <- "down.up"
fiafvc_volp$deg[abs(fiafvc_volp$log2FC) < log2(1.5)]              <- "ns"
fiafvc_volp <- fiafvc_volp[!fiafvc_volp$ens %in% DEG_variant, ]

fiafvc_volp$shared <- NA
fiafvc_volp$shared[fiafvc_volp$gene %in% c("Gm15501")] <- "Gm15501"
fiafvc_volp$shared[fiafvc_volp$gene %in% c("Htra1")]   <- "Htra1"
fiafvc_volp$shared[fiafvc_volp$gene %in% c("Cpne8")]   <- "Cpne8"
fiafvc_volp$shared[abs(fiafvc_volp$log2FC) < log2(1.5)] <- NA

fiafvc_volp_fr  <- which(fiafvc_volp$deg %in% c("up.up", "down.up"))
fiafvc_volp_beh <- setdiff(1:nrow(fiafvc_volp), fiafvc_volp_fr)
fiafvc_volp     <- fiafvc_volp[c(fiafvc_volp_beh, fiafvc_volp_fr), ]

ggplot(
  flfl_volp[, 2:6],
  aes(x = log2FC, y = -log10(p.adj), col = deg, label = shared)
) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "grey", "tomato", "darkgreen")) +
  geom_vline(xintercept = c(log2(2/3), log2(1.5)), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlim(-10, 18) +
  ylim(0, 18)

ggplot(
  fiafvc_volp[, 2:6],
  aes(x = log2FC, y = -log10(p.adj), col = deg, label = shared)
) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "grey", "tomato", "darkgreen")) +
  geom_vline(xintercept = c(log2(2/3), log2(1.5)), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlim(-10, 20) +
  ylim(0, 18)

#### Post-analysis counts ####

length(which(fiafvc_volp$log2FC < -0.6 & fiafvc_volp$p.adj < 0.05 & fiafvc_volp$deg == "sig"))
length(which(fiafvc_volp$log2FC > 0.6  & fiafvc_volp$p.adj < 0.05 & fiafvc_volp$deg == "sig"))
length(which(flfl_volp$log2FC > 0.6   & flfl_volp$p.adj < 0.05  & flfl_volp$deg == "sig"))
length(which(flfl_volp$log2FC < -0.6  & flfl_volp$p.adj < 0.05  & flfl_volp$deg == "sig"))

#### EnhancedVolcano plots ####

library(EnhancedVolcano)

EnhancedVolcano(
  flfl_volp,
  lab       = flfl_volp$gene,
  x         = "log2FC",
  y         = "p.adj",
  selectLab = c("Cpne8", "Htra1", "Gm15501"),
  pCutoff   = 0.05,
  FCcutoff  = log2(1.5),
  pointSize = 2.0,
  labSize   = 6.0,
  labCol    = "black",
  col       = c("grey", "grey", "grey", "red3"),
  labFace   = "bold",
  boxedLabels = TRUE,
  colAlpha    = 4/5,
  legendPosition = "right",
  legendLabSize  = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors  = "black"
)

down.up <- c("Cpne8", "Gm15501")
up.up   <- c("Htra1")

keyvals.colour <- ifelse(
  flfl_volp$deg == "ns", "grey",
  ifelse(
    flfl_volp$deg == "sig", "red",
    ifelse(flfl_volp$deg == "up.up", "blue", "darkgreen")
  )
)
names(keyvals.colour)[keyvals.colour == "grey"]      <- "n.s."
names(keyvals.colour)[keyvals.colour == "red3"]      <- "sig"
names(keyvals.colour)[keyvals.colour == "blue"]      <- "up.up"
names(keyvals.colour)[keyvals.colour == "darkgreen"] <- "down.up"

p_flfl <- EnhancedVolcano(
  flfl_volp,
  lab       = flfl_volp$gene,
  x         = "log2FC",
  y         = "p.adj",
  selectLab = c("Cpne8", "Htra1", "Gm15501"),
  xlab      = bquote(~Log[2]~"fold change"),
  title     = "flfl",
  pCutoff   = 0.05,
  FCcutoff  = log2(1.5),
  pointSize = 1.0,
  labSize   = 0.0,
  colCustom = keyvals.colour,
  colAlpha  = 1,
  legendPosition = "right",
  legendLabSize  = 15,
  legendIconSize = 5.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors  = "grey50",
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  border      = "full",
  borderWidth = 1.0,
  borderColour = "black"
)

#### GSEA input gene lists ####

flfl_genelist <- data.frame(
  gene  = flfl_volp$ens,
  value = flfl_volp$log2FC
)
flfl_genelist <- na.omit(flfl_genelist)
flfl_genelist <- flfl_genelist[order(flfl_genelist$value, decreasing = TRUE), ]
flfl_genelist <- flfl_genelist[!duplicated(flfl_genelist), ]
flfl_genelist1 <- flfl_genelist$value
names(flfl_genelist1) <- flfl_genelist$gene
write.csv(flfl_genelist, "genelist_flfl.csv")

fiafvc_genelist <- data.frame(
  gene  = fiafvc_volp$ens,
  value = fiafvc_volp$log2FC
)
fiafvc_genelist <- na.omit(fiafvc_genelist)
fiafvc_genelist <- fiafvc_genelist[order(fiafvc_genelist$value, decreasing = TRUE), ]
fiafvc_genelist <- fiafvc_genelist[!duplicated(fiafvc_genelist), ]
fiafvc_genelist1 <- fiafvc_genelist$value
names(fiafvc_genelist1) <- fiafvc_genelist$gene
write.csv(fiafvc_genelist, "genelist_fiafvc.csv")

##### Export for network analysis #####

flfl_volp1 <- data.frame(
  ens   = t.res.flfl$gene,
  gene  = t.res.flfl$symbol,
  log2FC = t.res.flfl$estimate,
  pval   = t.res.flfl$p.value,
  p.adj  = t.res.flfl$p.adjusted
)
write.csv(flfl_volp1, "genelist_flfl_annot.csv")

fiafvc_volp1 <- data.frame(
  ens   = t.res.fiafvc$gene,
  gene  = t.res.fiafvc$symbol,
  log2FC = t.res.fiafvc$estimate,
  pval   = t.res.fiafvc$p.value,
  p.adj  = t.res.fiafvc$p.adjusted
)
write.csv(fiafvc_volp1, "genelist_fiafvc_annot.csv")

#### msigdbr gene sets and GSEA ####

c5_mm <- msigdbr(species = "Mus musculus", category = "C5") %>%
  dplyr::select(gs_name, ensembl_gene)
gobp_mm <- msigdbr(
  species   = "Mus musculus",
  category  = "C5",
  subcategory = "GO:BP"
) %>%
  dplyr::select(gs_name, ensembl_gene)
c2_mm <- msigdbr(species = "Mus musculus", category = "C2") %>%
  dplyr::select(gs_name, ensembl_gene)

flfl_gsea <- GSEA(
  flfl_genelist1,
  TERM2GENE      = gobp_mm,
  pvalueCutoff   = 1,
  pAdjustMethod  = "fdr",
  seed           = 8888
)
write.csv(
  data.frame(
    ID   = flfl_gsea@result$ID,
    NES  = flfl_gsea@result$NES,
    p.adj = flfl_gsea@result$p.adjust,
    q.val = flfl_gsea@result$qvalue
  ),
  file = "flfl_gsea.csv"
)

fiafvc_gsea <- GSEA(
  fiafvc_genelist1,
  TERM2GENE      = gobp_mm,
  pvalueCutoff   = 1,
  pAdjustMethod  = "fdr",
  seed           = 8888
)
write.csv(
  data.frame(
    ID   = fiafvc_gsea@result$ID,
    NES  = fiafvc_gsea@result$NES,
    p.adj = fiafvc_gsea@result$p.adjust,
    q.val = fiafvc_gsea@result$qvalue
  ),
  file = "fiafvc_gsea.csv"
)

#### Dotplot combining selected GO terms for flfl and fiafvc ####

flfl_gsea_res <- flfl_gsea@result
flfl_gsea_filtered <- flfl_gsea[
  c(
    "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION",
    "GOBP_B_CELL_ACTIVATION",
    "GOBP_LEUKOCYTE_PROLIFERATION",
    "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM",
    "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
    "GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY",
    "GOBP_NEGATIVE_REGULATION_OF_CELL_CELL_ADHESION",
    "GOBP_MORPHOGENESIS_OF_AN_EPITHELIUM",
    "GOBP_COLLAGEN_METABOLIC_PROCESS",
    "GOBP_COMPLEMENT_ACTIVATION",
    "GOBP_CELLULAR_RESPONSE_TO_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_STIMULUS",
    "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION",
    "GOBP_POSITIVE_REGULATION_OF_BROWN_FAT_CELL_DIFFERENTIATION",
    "GOBP_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_ACUTE_PHASE_RESPONSE",
    "GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM",
    "GOBP_EPITHELIAL_STRUCTURE_MAINTENANCE",
    "GOBP_NEGATIVE_REGULATION_OF_PEPTIDASE_ACTIVITY"
  ),
  c("ID", "NES", "p.adjust")
]
flfl_gsea_filtered <- flfl_gsea_filtered[order(flfl_gsea_filtered$NES, decreasing = TRUE), ]
flfl_gsea_filtered$genotype <- "flfl"
flfl_gsea_filtered1 <- flfl_gsea_filtered[order(flfl_gsea_filtered$NES, decreasing = FALSE), ]

fiafvc_gsea_res <- fiafvc_gsea@result
fiafvc_gsea_filtered <- fiafvc_gsea[
  c(
    "GOBP_POSITIVE_REGULATION_OF_CELL_ACTIVATION",
    "GOBP_B_CELL_ACTIVATION",
    "GOBP_LEUKOCYTE_PROLIFERATION",
    "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM",
    "GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS",
    "GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY",
    "GOBP_NEGATIVE_REGULATION_OF_CELL_CELL_ADHESION",
    "GOBP_MORPHOGENESIS_OF_AN_EPITHELIUM",
    "GOBP_COLLAGEN_METABOLIC_PROCESS",
    "GOBP_COMPLEMENT_ACTIVATION",
    "GOBP_CELLULAR_RESPONSE_TO_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_STIMULUS",
    "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION",
    "GOBP_POSITIVE_REGULATION_OF_BROWN_FAT_CELL_DIFFERENTIATION",
    "GOBP_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_ACUTE_PHASE_RESPONSE",
    "GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM",
    "GOBP_EPITHELIAL_STRUCTURE_MAINTENANCE",
    "GOBP_NEGATIVE_REGULATION_OF_PEPTIDASE_ACTIVITY"
  ),
  c("ID", "NES", "p.adjust")
]
fiafvc_gsea_filtered <- fiafvc_gsea_filtered[order(fiafvc_gsea_filtered$NES, decreasing = TRUE), ]
fiafvc_gsea_filtered$genotype <- "fiafvc"

gsea_dotplot <- rbind(flfl_gsea_filtered, fiafvc_gsea_filtered)
gsea_dotplot$log10p <- -log10(gsea_dotplot$p.adjust)
gsea_dotplot$genotype <- factor(gsea_dotplot$genotype, levels = c("flfl", "fiafvc"))
gsea_dotplot$ID <- factor(gsea_dotplot$ID, levels = flfl_gsea_filtered1$ID)
gsea_dotplot$ID <- gsub("GOBP_", "", gsea_dotplot$ID)
gsea_dotplot$ID <- gsub("_", " ", gsea_dotplot$ID)

ID_levels <- gsub("GOBP_", "", flfl_gsea_filtered1$ID)
ID_levels <- gsub("_", " ", ID_levels)
gsea_dotplot$ID <- factor(gsea_dotplot$ID, levels = ID_levels)

library(ggthemes)
ggplot(
  gsea_dotplot,
  aes(x = genotype, y = ID, color = as.numeric(NES), size = as.numeric(log10p))
) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(colour = "NES", size = "-log10(fdr)") +
  theme_bw()

#### Additional GSEA plotting for flfl only ####

flfl_gsea <- GSEA(
  flfl_genelist1,
  TERM2GENE     = gobp_mm,
  pvalueCutoff  = 1,
  pAdjustMethod = "fdr",
  seed          = 8888
)

flfl_gsego <- gseGO(
  flfl_genelist1,
  ont            = "ALL",
  OrgDb          = org.Mm.eg.db,
  keyType        = "ENSEMBL",
  minGSSize      = 8,
  maxGSSize      = 500,
  pvalueCutoff   = 1,
  pAdjustMethod  = "fdr",
  seed           = 8888
)
write.csv(
  data.frame(
    ID   = flfl_gsego@result$Description,
    NES  = flfl_gsego@result$NES,
    p.adj = flfl_gsego@result$pvalue,
    q.val = flfl_gsego@result$qvalue
  ),
  file = "flfl_gsego.csv"
)

gseaplot2(
  flfl_gsego,
  geneSetID = c(
    291,  # negative cell-cell adhesion
    820,  # defense response to bacterium
    709,  # MyD88 TLR
    778,  # TLR
    28,   # FA metabolism
    38,
    62    # Regulation IL-6
  ),
  title = "2 weeks"
)

gseaplot2(
  flfl_gsego,
  geneSetID    = c(291, 820, 709, 778, 28, 38, 62),
  pvalue_table = TRUE,
  ES_geom      = "line"
)

############################################################
# End of script
############################################################