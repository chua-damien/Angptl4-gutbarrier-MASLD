############################################################
# Liver RNA-seq Analysis: Reversion + Inulin regimen
############################################################

#### Setup ####

## Core packages
library(DESeq2)
library(GEOquery)
library(gtools)
library(readr)
library(biobroom)
library(dplyr)
library(ggplot2)
library(ggforce)
library(pheatmap)
library(vsn)
library(sva)
library(tibble)
library(tidyverse)

## Enrichment & annotation
library(clusterProfiler)
library(AnnotationHub)
library(annotables)
library(org.Mm.eg.db)
library(msigdbr)

## ssGSEA
library(ssGSEA2)

set.seed(8888)

############################################################
#### 1. Load & Preprocess featureCounts Data ####
############################################################

temp_LPIR <- read.delim("./Datasets/fc_LPIR.txt")

# Rename all columns explicitly
colnames(temp_LPIR) <- c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length",
  # Batch 1: reversion & early LP samples
  "NAFL_R_1",   "NAFL_R_2",
  "LP_w8_r1",   "LP_w8_r2",
  "NASH_R_1",   "NASH_R_2",
  "Ctrl_w12_r1",
  # Batch 2: inulin samples
  "LPI_w12_1",  "LPI_w12_2",  "LPI_w12_3",
  "LP_w12_i1",  "LP_w12_i2",
  "LPI_w8_1",   "LPI_w8_2",
  # Batch 3: week 12
  "Ctrl_w12_1", "Ctrl_w12_2", "Ctrl_w12_3", "Ctrl_w12_4",
  "Ctrl_w12_5", "Ctrl_w12_6", "Ctrl_w12_7", "Ctrl_w12_8",
  "LP_w12_1",   "LP_w12_2",   "LP_w12_3",   "LP_w12_4",
  "LP_w12_5",   "LP_w12_6",   "LP_w12_7",   "LP_w12_8",
  # Batch 3: week 4
  "Ctrl_w4_1",  "Ctrl_w4_2",  "Ctrl_w4_3",  "Ctrl_w4_4",
  "Ctrl_w4_5",  "Ctrl_w4_6",  "Ctrl_w4_7",  "Ctrl_w4_8",
  "LP_w4_1",    "LP_w4_2",    "LP_w4_3",    "LP_w4_4",
  "LP_w4_5",    "LP_w4_6",    "LP_w4_7",    "LP_w4_8",
  # Batch 3: week 8
  "Ctrl_w8_1",  "Ctrl_w8_2",  "Ctrl_w8_3",  "Ctrl_w8_4",
  "Ctrl_w8_5",  "Ctrl_w8_6",  "Ctrl_w8_7",  "Ctrl_w8_8",
  "LP_w8_1",    "LP_w8_2",    "LP_w8_3",    "LP_w8_4",
  "LP_w8_5",    "LP_w8_6",    "LP_w8_7",    "LP_w8_8"
)

# Extract count matrix (gene ID + sample columns only)
LPIR_comb <- temp_LPIR[, c(1, 7:68)]

rownames(LPIR_comb) <- LPIR_comb[, 1]
LPIR_comb1          <- LPIR_comb[, 2:ncol(LPIR_comb)]

############################################################
#### 2. Build Sample Metadata ####
############################################################

LPIR_comb1_conds <- data.frame(
  SampleID = colnames(LPIR_comb1),
  batch = factor(c(rep(1, 7), rep(2, 7), rep(3, 48))),
  conditions = factor(
    c(
      rep("NAFL_R",    2), rep("w8_LP",    2),
      rep("NASH_R",    2), rep("w12_ctrl", 1),
      rep("w12_LPI",   3), rep("w12_LP",   2),
      rep("w8_LPI",    2),
      rep("w12_ctrl",  8), rep("w12_LP",   8),
      rep("w4_ctrl",   8), rep("w4_LP",    8),
      rep("w8_ctrl",   8), rep("w8_LP",    8)
    ),
    levels = c(
      "w4_ctrl",  "w4_LP",
      "w8_ctrl",  "w8_LP",
      "w12_ctrl", "w12_LP",
      "NAFL_R",   "NASH_R",
      "w8_LPI",   "w12_LPI"
    )
  ),
  row.names = NULL
)
rownames(LPIR_comb1_conds) <- LPIR_comb1_conds$SampleID


############################################################
#### 3. DESeq2 (ComBat-seq Batch Correction) ####
############################################################

covar_mat           <- as.matrix(LPIR_comb1_conds[, c("batch", "conditions")])
rownames(covar_mat) <- LPIR_comb1_conds$SampleID

LPIR_cbs <- ComBat_seq(
  as.matrix(LPIR_comb1),
  batch = covar_mat[, "batch"],
  group = covar_mat[, "conditions"]
)

model.matrix(~ conditions, LPIR_comb1_conds)

dds_LPIR_cbs <- DESeqDataSetFromMatrix(
  countData = as.matrix(LPIR_cbs),
  colData   = as.data.frame(LPIR_comb1_conds),
  design    = ~ batch + conditions
)
dds_LPIR_cbs <- DESeq(dds_LPIR_cbs)
resultsNames(dds_LPIR_cbs)

res.LPIR_cbs   <- results(dds_LPIR_cbs)
t.res.LPIR_cbs <- tidy.DESeqResults(res.LPIR_cbs)
t.res.LPIR_cbs <- arrange(t.res.LPIR_cbs, p.adjusted)

#### 4a. VST & PCA (with BC) ####

vds_LPIR_cbs <- vst(dds_LPIR_cbs, blind = FALSE)
vsn::meanSdPlot(assay(vds_LPIR_cbs))

plotPCA(vds_LPIR_cbs, intgroup = "conditions") +
  ggforce::geom_mark_ellipse(aes(color = conditions))

pca_LPIR_cbs <- plotPCA(vds_LPIR_cbs, intgroup = "conditions", returnData = TRUE)
percentVar   <- round(100 * attr(pca_LPIR_cbs, "percentVar"))

ggplot(pca_LPIR_cbs, aes(PC1, PC2, color = conditions, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  geom_text(hjust = 0, vjust = 0) +
  ggforce::geom_mark_ellipse(aes(fill = conditions, color = conditions))

# Shared column annotation for heatmaps (defined once here for reuse below)
df_col <- data.frame(conditions = colData(dds_LPIR_cbs)[, "conditions"],
                     row.names  = colData(dds_LPIR_cbs)[, "SampleID"])

annot_col <- list(conditions = c(
  w4_ctrl  = "black",  w12_ctrl = "black",
  w8_LP    = "red",    w12_LP   = "red",
  NAFL_R   = "lightblue",
  w8_LPI   = "green",  w12_LPI  = "green"
))

############################################################
#### 4. VST Export for ssGSEA ####
############################################################
saveRDS(vds_LPIR_cbs, "LPIR_vst.rds")

vds_LPIR_cbs.mat <- assay(vds_LPIR_cbs)
cat("VST matrix dimensions:", dim(vds_LPIR_cbs.mat), "\n")

vds_LPIR_cbs.mat <- tibble::rownames_to_column(as.data.frame(vds_LPIR_cbs.mat), "gene")
vds_LPIR_cbs.mat <- inner_join(vds_LPIR_cbs.mat, grcm38, by = c("gene" = "ensgene"))

# Replace Ensembl IDs with UPPERCASE gene symbols (GSEA/GCT format)
vds_LPIR_cbs.ex      <- vds_LPIR_cbs.mat
vds_LPIR_cbs.ex[, 1] <- toupper(vds_LPIR_cbs.mat$symbol)
vds_LPIR_cbs.ex      <- vds_LPIR_cbs.ex[, 1:63]  # gene symbol + 62 sample columns

write.table(
  vds_LPIR_cbs.ex,
  "LPIR_vst.txt",
  row.names = FALSE,
  col.names = TRUE,
  sep       = "\t"
)


############################################################
#### 5. ssGSEA – NAFLD Signature v2 (Human Gene Set Merged) ####
############################################################

vds_LPIR_cbs.ex.2 <- distinct(vds_LPIR_cbs.ex, gene, .keep_all = TRUE)

vds_meta           <- read.csv("./ssGSEA-MASLDsig-v2/hs.pkg/vsd_all2-damien.csv")
colnames(vds_meta)[1] <- "gene"

full_gene <- intersect(vds_LPIR_cbs.ex.2$gene, vds_meta$gene)
vds_all   <- vds_LPIR_cbs.ex.2[vds_LPIR_cbs.ex.2$gene %in% full_gene, ]
vds_all   <- merge(vds_all, vds_meta, by = "gene")

write.table(vds_all, "LPIR_vds_all.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

res_ssgsea_nafld_v2 <- run_ssGSEA2(
  "./LPIR_vds_all-reformat.gct",
  output.prefix      = "LPIR.vsd",
  gene.set.databases = "./ssGSEA-MASLDsig-v2/hs.pkg/NAFLDgeneset.gmt",
  output.directory   = "./ssGSEA-MASLDsig-v2/",
  sample.norm.type   = "none",
  weight             = 0.75,
  correl.type        = "rank",
  statistic          = "area.under.RES",
  output.score.type  = "NES",
  nperm              = 1000,
  min.overlap        = 5,
  extended.output    = TRUE,
  global.fdr         = FALSE,
  log.file           = "./run.log"
)

############################################################
#### 6. ssGSEA – GO:BP Selected Genesets ####
############################################################

res_ssgsea_gobp <- run_ssGSEA2(
  "./LPIR_vst-reformat.gct",
  output.prefix      = "LPIR.vsd.gobp_selected",
  gene.set.databases = "./ssGSEA-genesets/genesets/GOBP_selected.gmt",
  output.directory   = "./ssGSEA-genesets/",
  sample.norm.type   = "none",
  weight             = 0.75,
  correl.type        = "rank",
  statistic          = "area.under.RES",
  output.score.type  = "NES",
  nperm              = 1000,
  min.overlap        = 5,
  extended.output    = TRUE,
  global.fdr         = FALSE,
  log.file           = "./run.log",
  spare.cores        = 16,
  par                = FALSE
)

############################################################
#### 7. ssGSEA Output Heatmaps ####
############################################################

ssgsea_res <- read.csv("./ssGSEA-genesets/LPIR_gobp-scores.csv",
                        header = TRUE, row.names = 1)
rownames(ssgsea_res) <- gsub("GOBP_", "", rownames(ssgsea_res))

# Shared pathway terms for focused heatmaps
SELECTED_TERMS <- c(
  "FATTY_ACID_METABOLIC_PROCESS",
  "FATTY_ACID_BETA_OXIDATION",
  "ORGANIC_HYDROXY_COMPOUND_BIOSYNTHETIC_PROCESS",
  "ORGANIC_HYDROXY_COMPOUND_CATABOLIC_PROCESS",
  "LIPID_BIOSYNTHETIC_PROCESS",
  "REGULATION_OF_LIPID_LOCALIZATION",
  "LIPID_METABOLIC_PROCESS",
  "ACUTE_INFLAMMATORY_RESPONSE",
  "ACUTE_PHASE_RESPONSE",
  "TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "CHRONIC_INFLAMMATORY_RESPONSE",
  "ENDOTHELIAL_CELL_PROLIFERATION",
  "ENDOTHELIAL_CELL_MATRIX_ADHESION",
  "EXTRACELLULAR_MATRIX_DISASSEMBLY",
  "COLLAGEN_CATABOLIC_PROCESS",
  "COLLAGEN_METABOLIC_PROCESS",
  "COLLAGEN_FIBRIL_ORGANIZATION"
)

MASL_R_SAMPLES <- c("Ctrl_w4_1", "Ctrl_w4_2",
                     "LP_w8_r1",  "LP_w8_r2",
                     "NAFL_R_1",  "NAFL_R_2")

INULIN_SAMPLES <- c(
  "Ctrl_w12_3", "Ctrl_w12_4", "Ctrl_w12_6", "Ctrl_w12_7", "Ctrl_w12_8",
  "LP_w12_i1",  "LP_w12_i2",  "LP_w12_6",   "LP_w12_7",   "LP_w12_8",
  "LPI_w8_1",   "LPI_w8_2",   "LPI_w12_1",  "LPI_w12_2",  "LPI_w12_3"
)

HEATMAP_COLORS <- colorRampPalette(c("blue", "white", "red"))(75)

#### 7a. MASL_R: All pathways (clustered) ####

pheatmap(
  ssgsea_res[, MASL_R_SAMPLES],
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  annotation_col    = df_col,
  annotation_colors = annot_col,
  main              = "ssGSEA: MASL_R (all pathways)"
)

#### 7b. Inulin: All pathways (clustered) ####

pheatmap(
  ssgsea_res[, INULIN_SAMPLES],
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  annotation_col    = df_col,
  annotation_colors = annot_col,
  main              = "ssGSEA: Inulin (all pathways)"
)


############################################################
# End of script
############################################################