############################################################
# Fecal 16S rRNA-seq Analysis: Diversity + Gram status pieline + Composition + Metabolic pathways
# Sections:
#   1.  Build Phyloseq Object
#   2.  ConQuR Batch Correction
#   3.  Alpha Diversity
#   4.  Beta Diversity (Ordination, PERMANOVA)
#   5.  Top Taxa Bar Plots
#   6.  Dysbiosis Score
#   7.  Individual Bacteria Plots (Absolute & Relative)
#   8.  Gram Identity Database Building
#   9.  Tax Annotation (Gram Mapping)
#   10. Gram Ratio Calculation
#   11. Top Gram-Stratified Bar Plots 
#   12. Metabolic Pathway analysis
############################################################

#### Setup ####

## Core microbiome packages
library(phyloseq)
library(Biostrings)
library(dada2)
library(ape)
library(vegan)
library(microViz)
library(microbiome)
library(microshades)
library(VennDiagram)

## Batch correction
library(ConQuR)
library(foreach)

## Gram identity & annotation
library(epidm)
library(stringr)

## Diversity & dysbiosis
library(dysbiosisR)

## Plotting & wrangling
library(ggplot2)
library(ggforce)
library(ggpubr)
library(dplyr)
library(data.table)
library(matrixStats)

theme_set(theme_bw())

#### Helper Functions ####

# Remove specific taxa from a phyloseq object
pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

# Helper: build epsortdb subset by GramStain category
filter_gram <- function(db, stain_value) {
  cols_keep <- c("TaxID", "Organism", "Genus", "Species", "Phylum", "Class", "GramStain")
  out <- db[db$GramStain %in% stain_value, which(colnames(db) %in% cols_keep)]
  return(unique(out))
}

############################################################
#### 1. Build Phyloseq Object ####
############################################################

samples.out <- rownames(seqtab.nochim)
sample.id   <- gsub("-", "", substr(samples.out, 0, 5))

samdf <- data.frame(
  Diet     = factor(c("LIDPAD","LIDPAD","Ctrl","Ctrl","LIDPAD","LIDPAD","LIDPAD","LIDPAD","Ctrl","Ctrl","Ctrl","Ctrl"),
                    levels = c("Ctrl", "LIDPAD")),
  Genotype = factor(c("Angptl4.vc","Angptl4.flfl","Angptl4.vc","Angptl4.flfl",
                       "Angptl4.flfl","Angptl4.vc","Angptl4.flfl","Angptl4.vc",
                       "Angptl4.vc","Angptl4.vc","Angptl4.flfl","Angptl4.flfl")),
  Cage     = factor(c(1,1,2,2,3,3,3,3,4,4,4,4)),
  SampleID = sample.id
)
rownames(samdf) <- samples.out

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa)
)

# Attach DNA sequences and rename ASVs
dna          <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna)   <- taxa_names(ps)
ps           <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Remove known contaminant ASV
tax_table(ps)[which(rownames(tax_table(ps)) %in% "ASV17631"), ]
ps <- pop_taxa(ps, "ASV17631")

# Add random tree for UniFrac
random_tree <- ape::rtree(ntaxa(ps), rooted = TRUE, tip.label = taxa_names(ps))
ps          <- merge_phyloseq(ps, random_tree)

saveRDS(ps, "LG_vilcre_ps.rds")

#### Load saved phyloseq ####
ps <- readRDS("LG_vilcre_ps.rds")

############################################################
#### 2. ConQuR Batch Correction ####
############################################################

taxa_raw <- otu_table(ps)
batchid  <- sample_data(ps)$Cage
covar    <- sample_data(ps)[, c("Diet", "Genotype")]

taxa_corrected1 <- ConQuR(tax_tab = taxa_raw, batchid = batchid,
                           covariates = covar, batch_ref = "4")

options(warn = -1)
taxa_corrected2 <- ConQuR(tax_tab = taxa_raw, batchid = batchid,
                           covariates = covar, batch_ref = "4",
                           logistic_lasso = TRUE, quantile_type = "lasso", interplt = TRUE)

# Visualise batch effect correction
par(mfrow = c(2, 3))
Plot_PCoA(TAX = taxa_raw,         factor = batchid, main = "Before Correction, Bray-Curtis")
Plot_PCoA(TAX = taxa_corrected1,  factor = batchid, main = "ConQuR (Default), Bray-Curtis")
Plot_PCoA(TAX = taxa_corrected2,  factor = batchid, main = "ConQuR (Penalized), Bray-Curtis")
Plot_PCoA(TAX = taxa_raw,         factor = batchid, dissimilarity = "Aitch", main = "Before Correction, Aitchison")
Plot_PCoA(TAX = taxa_corrected1,  factor = batchid, dissimilarity = "Aitch", main = "ConQuR (Default), Aitchison")
Plot_PCoA(TAX = taxa_corrected2,  factor = batchid, dissimilarity = "Aitch", main = "ConQuR (Penalized), Aitchison")
par(mfrow = c(1, 1))

# Apply corrected OTU table (use penalized version)
ps_cqn <- ps
otu_table(ps_cqn) <- otu_table(taxa_corrected2, taxa_are_rows = FALSE)


############################################################
#### 3. Alpha Diversity ####
############################################################

ALPHA_MEASURES <- c("Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed")
DIET_COLORS    <- c("black", "red")

# All samples: overview across genotypes
plot_richness(ps_cqn_omit, x = "Diet",
              measures = c("Observed", "Shannon", "Simpson", "Chao1"),
              color = "Genotype") +
  geom_point(size = 3)

# flfl only: full metric panel
plot_richness(ps_cqn_flfl, x = "Diet", measures = ALPHA_MEASURES, color = "Diet") +
  geom_boxplot(alpha = 0.1, width = 0.75) +
  geom_point(size = 2, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = DIET_COLORS)

# Statistical tests (Kruskal + Wilcoxon, one-sided: LIDPAD < Ctrl)
erich_omit <- data.frame(sample_data(ps_cqn_flfl),
                          estimate_richness(ps_cqn_flfl, measures = ALPHA_MEASURES))

for (metric in ALPHA_MEASURES) {
  cat("\n---", metric, "---\n")
  cat("Kruskal-Wallis:\n")
  print(kruskal.test(reformulate("Diet", response = metric), data = erich_omit))
  cat("Wilcoxon (LIDPAD < Ctrl):\n")
  print(wilcox.test(
    erich_omit[[metric]][erich_omit$Diet == "LIDPAD"],
    erich_omit[[metric]][erich_omit$Diet == "Ctrl"],
    alternative = "less", paired = FALSE, exact = FALSE, correct = FALSE
  ))
}

############################################################
#### 4. Beta Diversity: Ordination & PERMANOVA ####
############################################################

# Proportional transform
ps.prop      <- transform_sample_counts(ps_cqn_omit, function(otu) otu / sum(otu))
ps.prop.flfl <- transform_sample_counts(ps_cqn_flfl, function(otu) otu / sum(otu))

# NMDS Bray-Curtis (all samples)
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")

plot_ordination(ps.prop, ord.nmds.bray, color = "Diet", title = "Bray NMDS") +
  facet_grid(~ Genotype) +
  stat_ellipse(type = "norm", linetype = 1) +
  expand_limits(x = c(-0.4, 0.4), y = c(-0.45, 0.5)) +
  ggforce::geom_mark_ellipse(aes(color = Diet)) +
  scale_color_manual(values = DIET_COLORS)

# PCoA Bray-Curtis (flfl only)
ord.pcoa.bray.flfl <- ordinate(ps.prop.flfl, method = "PCoA", distance = "bray")

plot_ordination(ps.prop.flfl, ord.pcoa.bray.flfl, color = "Diet", title = "Bray PCoA (flfl)") +
  ggforce::geom_mark_ellipse(aes(color = Diet)) +
  expand_limits(x = c(-0.5, 0.6), y = c(-0.35, 0.4)) +
  scale_color_manual(values = DIET_COLORS)

# microViz interactive explorer (run interactively)
# ord_explore(ps_cqn_flfl)

# PERMANOVA across multiple distances
D_BC  <- phyloseq::distance(ps.prop.flfl, "bray")
D_JC  <- phyloseq::distance(ps.prop.flfl, "jaccard")
D_UF  <- phyloseq::distance(ps.prop.flfl, "unifrac")
D_wUF <- phyloseq::distance(ps.prop.flfl, "wunifrac")

dist_list <- list(
  "Bray Curtis"     = D_BC,
  "Jaccard"         = D_JC,
  "UniFrac"         = D_UF,
  "Weighted UniFrac" = D_wUF
)

meta_flfl <- as(sample_data(ps_cqn_flfl), "data.frame")

for (i in seq_along(dist_list)) {
  res <- vegan::adonis2(dist_list[[i]] ~ Diet, data = meta_flfl)
  cat(names(dist_list[i]), "\n")
  print(res)
  cat("\n", strrep("-", 70), "\n")
}

############################################################
#### 5. Top Taxa Bar Plots ####
############################################################

plot_top_taxa <- function(ps_obj, rank_level, n_taxa, title_prefix = "") {
  ps_trf <- tax_fix(ps_obj)
  ps_trf <- tax_transform(ps_trf, trans = "compositional", rank = rank_level)

  top_agg    <- aggregate_taxa(ps_obj, level = rank_level)
  top_names  <- top_taxa(top_agg, n = n_taxa + 1)
  top_names  <- top_names[2:(n_taxa + 1)]  # exclude "Unknown"
  top_pruned <- prune_taxa(top_names, top_agg)
  top_pruned <- transform_sample_counts(top_pruned, function(ASV) ASV / sum(ASV))

  plot_bar(top_pruned, x = "SampleID", fill = rank_level) +
    facet_wrap(~ Diet + Genotype, scales = "free_x", ncol = 5) +
    scale_fill_manual(values = c(
      microshades_palette("micro_cvd_gray"),
      microshades_palette("micro_cvd_blue"),
      microshades_palette("micro_cvd_orange"),
      microshades_palette("micro_brown")
    )) +
    ggtitle(paste0(title_prefix, " Top ", n_taxa, " (", rank_level, ")"))
}

# All samples - Class level
plot_top_taxa(ps_cqn, "Class", 20, "All samples")

# flfl only - Phylum level
ps_cqn_omit1_phy <- subset_samples(ps_cqn, SampleID %in% SAMPLES_OMIT)
ps_cqn_flfl1_phy <- subset_samples(ps_cqn_omit1_phy, Genotype == "Angptl4.flfl")
plot_top_taxa(ps_cqn_flfl1_phy, "Phylum", 20, "flfl")
plot_top_taxa(ps_cqn_flfl1_phy, "Phylum", 50, "flfl")

############################################################
#### 6. Dysbiosis Score ####
############################################################
# Ref: https://microsud.github.io/dysbiosisR/articles/Introduction.html

ps_trf_dys <- tax_fix(subset_samples(ps_cqn_omit, Genotype == "Angptl4.flfl"))
ps_trf_dys <- tax_transform(ps_trf_dys, trans = "compositional", rank = "Genus")

dist.mat    <- phyloseq::distance(ps_trf_dys, "bray")
ref.samples <- sample_names(subset_samples(ps_trf_dys, Diet == "Ctrl"))

dysbiosis_1 <- dysbiosisMedianCLV(ps_trf_dys,
                                   dist_mat          = dist.mat,
                                   reference_samples = ref.samples)

dysbiosis_thres    <- quantile(subset(dysbiosis_1, Diet == "LIDPAD")$score, 0.9)
normobiosis_thres  <- quantile(subset(dysbiosis_1, Diet == "LIDPAD")$score, 0.1)

dysbiosis_1 <- dysbiosis_1 |>
  mutate(isDysbiostic = ifelse(score >= dysbiosis_thres, TRUE, FALSE))

ggplot(dysbiosis_1, aes(x = Diet, y = score)) +
  geom_boxplot(aes(middle = mean(score)), alpha = 0.1, width = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  ggtitle("DysbiosisR Score")

kruskal.test(score ~ Diet, data = dysbiosis_1)
t.test(score ~ Diet, data = dysbiosis_1, alternative = "two.sided", var.equal = TRUE)

############################################################
#### 7. Individual Bacteria Plots ####
############################################################

GENERA_OF_INTEREST <- c("Romboutsia", "Roseburia")

prepare_genus_df <- function(ps_obj, genera, relative = FALSE) {
  ps_trf <- subset_taxa(ps_obj, Kingdom != "unassigned")
  ps_trf <- tax_fix(ps_trf, unknowns = c("Incertae Sedis"), min_length = 2, suffix_rank = "Genus")
  if (relative) {
    ps_trf <- transform_sample_counts(ps_trf, function(ASV) ASV / sum(ASV))
  }
  ps_sel  <- subset_taxa(ps_trf, Genus %in% genera)
  ps_sel1 <- tax_glom(ps_sel, taxrank = "Genus")

  df       <- as.data.frame(otu_table(ps_sel1))
  df$SampleID <- rownames(df)
  sdf      <- as.data.frame(sample_data(ps_trf))
  sdf$SampleID <- rownames(sdf)
  df       <- dplyr::inner_join(df, sdf, by = "SampleID")
  df$Group <- paste(df$Diet, df$Genotype, sep = "_")

  # Identify ASV-to-genus mapping from tax table
  tt <- as.data.frame(tax_table(ps_sel1))
  asv_map <- setNames(tt$Genus, rownames(tt))

  df <- data.table::melt(
    df,
    id.vars      = c("SampleID", "Diet", "Genotype", "Cage", "Group"),
    variable.name = "ASV",
    direction    = "long",
    variable.factor = FALSE
  )
  df$ASV    <- as.character(df$ASV)
  df$Genus  <- asv_map[df$ASV]
  return(df)
}

# Absolute abundance
df_abs <- prepare_genus_df(ps_cqn, GENERA_OF_INTEREST, relative = FALSE)

ggplot(df_abs[df_abs$Genotype == "Angptl4.flfl", ], aes(x = Diet, y = value, color = Diet)) +
  geom_violin(outliers = FALSE, linewidth = 1) +
  geom_point(size = 1, position = position_jitter(w = 0.05)) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "black") +
  facet_wrap(~ Genus) +
  scale_color_manual(values = DIET_COLORS) +
  ggtitle("Absolute Abundance") +
  guides(color = "none")

# Relative abundance
df_rel <- prepare_genus_df(ps_cqn, GENERA_OF_INTEREST, relative = TRUE)

ggplot(df_rel[df_rel$Genotype == "Angptl4.flfl", ], aes(x = Diet, y = value)) +
  geom_boxplot(aes(middle = mean(value)), alpha = 0.1, width = 0.75) +
  geom_jitter(width = 0.05) +
  facet_wrap(~ Genus, scales = "free") +
  ggtitle("Relative Abundance") + labs(y = "Relative abundance")

# Statistical tests per genus
for (g in GENERA_OF_INTEREST) {
  cat("\n---", g, "---\n")
  sub <- df_rel[df_rel$Genus == g, ]
  print(kruskal.test(value ~ Diet, data = sub))
  print(wilcox.test(value ~ Diet, data = sub))
}

############################################################
#### 8. Gram Identity Database Building ####
############################################################
# Sources:
#   (1) PSORTdb v4: https://db.psort.org/
#   (2) epidm R package: genus_gram_stain

##### 8a. PSORTdb (epsortdb) #####

epsortdb  <- read.delim("../Databases/Experimental-PSORTdb-v4.00.tsv")

epsortdb1 <- epsortdb
epsortdb1[, c("Genus", "Species")] <- stringr::str_split_fixed(epsortdb1$Organism, " ", 2)

epsortdb1.trunc <- unique(epsortdb1[, which(colnames(epsortdb1) %in%
                                              c("TaxID", "Organism", "Genus", "Species", "Phylum", "Class", "GramStain"))])

# Subset by Gram stain category
epsort.gpos  <- filter_gram(epsortdb1, "Gram positive")
epsort.gneg  <- filter_gram(epsortdb1, "Gram negative")
epsort.gposv <- filter_gram(epsortdb1, "Gram positive with outer membrane")
epsort.gnegv <- filter_gram(epsortdb1, "Gram negative without outer membrane")  # Mycoplasma only

# Phylum Venn diagram (epsortdb)
VennDiagram::venn.diagram(
  x              = list(unique(epsort.gpos$Phylum), unique(epsort.gneg$Phylum)),
  category.names = c("gpos.phy", "gneg.phy"),
  filename       = "epsort.phy.png"
)

# Manual curation: correct mis-assigned genera
epsortdb1[epsortdb1$Genus %in% c("Staphylococcus", "Thermoanaerobacter",
                                    "Streptomyces", "Mycobacterium", "Corynebacterium"),
          "GramStain"] <- "Gram positive"
epsortdb1.trunc[epsortdb1.trunc$Genus %in% c("Staphylococcus", "Thermoanaerobacter",
                                               "Streptomyces"), "GramStain"] <- "Gram positive"

# Refresh after curation
epsort.gpos <- filter_gram(epsortdb1, "Gram positive")
epsort.gneg <- filter_gram(epsortdb1, "Gram negative")

##### 8b. epidm database #####

epidm.db <- genus_gram_stain
epidm.db$organism_genus <- stringr::str_to_title(epidm.db$organism_genus, locale = "en")
epidm.db$gram_stain[epidm.db$gram_stain == "POSITIVE"] <- "Gram positive"
epidm.db$gram_stain[epidm.db$gram_stain == "NEGATIVE"] <- "Gram negative"
colnames(epidm.db)[1:2] <- c("Genus", "GramStain")

epidm.gpos <- unique(epidm.db[epidm.db$GramStain == "Gram positive", ])
epidm.gneg <- unique(epidm.db[epidm.db$GramStain == "Gram negative", ])

##### 8c. Aggregated DB (agg.db) #####

# Genus-level Venn diagram comparing both databases
VennDiagram::venn.diagram(
  x = list(unique(epsort.gpos$Genus), unique(epsort.gneg$Genus),
            unique(epidm.gpos$Genus), unique(epidm.gneg$Genus)),
  category.names = c("db1.gpos.genus", "db1.gneg.genus",
                     "db2.gpos.genus", "db2.gneg.genus"),
  filename = "agg.db.png"
)

# Check conflicting genera across databases
cat("Conflicts in epsortdb:\n")
print(intersect(unique(epsort.gpos$Genus), unique(epsort.gneg$Genus)))
cat("epidm.gneg in epsort.gpos:\n")
print(intersect(unique(epidm.gneg$Genus), unique(epsort.gpos$Genus)))
cat("epidm.gpos in epsort.gneg:\n")
print(intersect(unique(epidm.gpos$Genus), unique(epsort.gneg$Genus)))

# Merge databases
agg.db <- dplyr::bind_rows(epsortdb1, epidm.db[, 1:2])

# Manual curation of known conflicts
agg.db[agg.db$Genus %in% c("Azospirillum", "Acidaminococcus", "Mycoplasma",
                               "Spiroplasma", "Parabacteroides"), "GramStain"] <- "Gram negative"
agg.db[agg.db$Phylum %in% "Firmicutes", "GramStain"] <- "Gram positive"

##### 8d. Genus-level reference DB (agg.db.genus) #####

agg.db.genus <- unique(agg.db[, c("Genus", "GramStain")])

# Remaining manual overrides
agg.db.genus[agg.db.genus$Genus %in% c("Candidatus", "UBA1819",
                                          "Oscillibacter",
                                          "Lachnospiraceae NK4A136 group"), "GramStain"] <- "Gram negative"
agg.db.genus[agg.db.genus$Genus %in% c("Deinococcus", "Staphylococcus",
                                          "Thermoanaerobacter", "Streptomyces",
                                          "Mycobacterium", "Corynebacterium"), "GramStain"] <- "Gram positive"

# Remove plasmid entries
agg.db.genus <- agg.db.genus[agg.db.genus$Genus != "Plasmid", ]
agg.db.genus <- unique(agg.db.genus)

cat("Gram status distribution (genus-level):\n")
print(table(agg.db.genus$GramStain))

##### 8e. Phylum-level reference DB (agg.db.phy) #####

agg.db.phy <- unique(agg.db[, c("Phylum", "GramStain")])
agg.db.phy[agg.db.phy$Phylum %in% "Deinococcus-Thermus", "GramStain"] <- "Gram negative"
agg.db.phy[agg.db.phy$Phylum %in% "Firmicutes",          "GramStain"] <- "Gram positive"
agg.db.phy <- agg.db.phy[!agg.db.phy$Phylum %in% c("", NA), ]
agg.db.phy <- unique(agg.db.phy)

# Extend with SILVA-style phylum names used in mouse gut 16S
agg.dbb.phy <- rbind(agg.db.phy, data.frame(
  Phylum    = c("Deferribacterota", "Bacteroidota", "Actinobacteriota",
                "Desulfobacterota", "Bacteria Kingdom", "Patescibacteria",
                "Campylobacterota", "Verrucomicrobiota", "Acidobacteriota"),
  GramStain = c("Gram negative", "Gram negative", "Gram positive",
                "Gram negative", "Unknown",        "Gram positive",
                "Gram negative", "Gram negative",  "Gram negative")
))


############################################################
#### 9. Tax Annotation (Gram Mapping) ####
############################################################

# Working phyloseq: flfl only, genus-level log transform
rank <- "Genus"

ps_cqn_omit1 <- subset_samples(ps_cqn, SampleID %in% SAMPLES_OMIT)
ps_cqn_flfl1 <- subset_samples(ps_cqn_omit1, Genotype == "Angptl4.flfl")
ps_trf        <- tax_fix(ps_cqn_flfl1, min_length = 2)
ps_trf        <- tax_transform(ps_trf, trans = "log", zero_replace = "halfmin", rank = rank)

cat("Coverage check — genera in ps_trf vs agg.db.genus:\n")
cat("  ps_trf genera:      ", length(colnames(otu_table(ps_trf))), "\n")
cat("  matched in agg.db:  ", length(intersect(colnames(otu_table(ps_trf)), agg.db.genus$Genus)), "\n")

## Genus-level tax annotation
tax_annot <- dplyr::inner_join(
  y  = agg.dbb.phy,
  x  = as.data.frame(tax_table(ps_trf)),
  by = "Phylum"
)
# Refine with genus-level assignments (genus overrides phylum)
tax_annot <- dplyr::left_join(y = agg.db.genus,      x = tax_annot, by = "Genus")
tax_annot <- dplyr::left_join(y = agg.db.genus.flag, x = tax_annot, by = "Genus", multiple = "first")

tax_annot.1 <- unique(tax_annot)
# Prefer genus-level GramStain (GramStain.y) over phylum-level (GramStain.x) where available
tax_annot.1[!is.na(tax_annot.1$GramStain.y), "GramStain.x"] <-
  tax_annot.1[!is.na(tax_annot.1$GramStain.y), "GramStain.y"]

tax_annot.1 <- tax_annot.1[, c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "GramStain.x", "Flagellin")]
colnames(tax_annot.1)[7] <- "GramStain"
rownames(tax_annot.1)    <- tax_annot.1$Genus

## ASV-level tax annotation
ps_trf_asv        <- tax_fix(ps_cqn_flfl1, min_length = 2)
tax_annot.asv     <- as.data.frame(tax_table(ps_trf_asv))
tax_annot.asv$ASV <- rownames(tax_annot.asv)

tax_annot.asv <- dplyr::left_join(y = agg.dbb.phy,      x = tax_annot.asv, by = "Phylum")
tax_annot.asv <- dplyr::left_join(y = agg.db.genus,      x = tax_annot.asv, by = "Genus")
tax_annot.asv <- dplyr::left_join(y = agg.db.genus.flag, x = tax_annot.asv, by = "Genus", multiple = "first")

tax_annot.asv.1 <- unique(tax_annot.asv)
tax_annot.asv.1[!is.na(tax_annot.asv.1$GramStain.y), "GramStain.x"] <-
  tax_annot.asv.1[!is.na(tax_annot.asv.1$GramStain.y), "GramStain.y"]

tax_annot.asv.1 <- tax_annot.asv.1[, c("Kingdom", "Phylum", "Class", "Order",
                                          "Family", "Genus", "ASV", "GramStain.x", "Flagellin")]
colnames(tax_annot.asv.1)[8] <- "GramStain"

############################################################
#### 10. Gram Ratio Calculation ####
############################################################

otu.count.asv <- otu_table(ps_trf_asv)

neg_asvs <- tax_annot.asv.1[tax_annot.asv.1$GramStain == "Gram negative", "ASV"]
pos_asvs <- tax_annot.asv.1[tax_annot.asv.1$GramStain == "Gram positive", "ASV"]

df_gram_ratio <- data.frame(
  gram.neg.counts = rowSums(otu.count.asv[, colnames(otu.count.asv) %in% neg_asvs]),
  gram.pos.counts = rowSums(otu.count.asv[, colnames(otu.count.asv) %in% pos_asvs])
)

df_gram_ratio$neg.pos   <- df_gram_ratio$gram.neg.counts / df_gram_ratio$gram.pos.counts
df_gram_ratio$total.cts <- df_gram_ratio$gram.neg.counts + df_gram_ratio$gram.pos.counts
df_gram_ratio$neg.comp  <- df_gram_ratio$gram.neg.counts / df_gram_ratio$total.cts
df_gram_ratio$pos.comp  <- df_gram_ratio$gram.pos.counts / df_gram_ratio$total.cts
df_gram_ratio           <- cbind(as.data.frame(sample_data(ps_trf_asv)), df_gram_ratio)

# Gram neg:pos ratio box-dot plot
ggplot(df_gram_ratio, aes(x = Diet, y = neg.pos)) +
  geom_boxplot(aes(middle = mean(neg.pos)), alpha = 0.1, width = 0.75) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  scale_fill_manual(values = DIET_COLORS) +
  ggtitle("Ratio of Gram Identity") + labs(y = "Gram neg : Gram pos")

# Stacked composition plot
df_gram_ratio.long <- data.table::melt(
  df_gram_ratio[, c("Diet", "SampleID", "neg.comp", "pos.comp")],
  id.vars       = c("SampleID", "Diet"),
  variable.name = "GramStain",
  direction     = "long"
)

ggplot(df_gram_ratio.long, aes(x = SampleID, y = value, fill = GramStain)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~ Diet, space = "free_x", scales = "free")

# Statistical tests
cat("t-test (Gram neg:pos ratio):\n")
print(t.test(
  df_gram_ratio$neg.pos[df_gram_ratio$Diet == "Ctrl"],
  df_gram_ratio$neg.pos[df_gram_ratio$Diet == "LIDPAD"],
  alternative = "two.sided", var.equal = TRUE
))

cat("Wilcoxon test:\n")
print(wilcox.test(
  df_gram_ratio$neg.pos[df_gram_ratio$Diet == "Ctrl"],
  df_gram_ratio$neg.pos[df_gram_ratio$Diet == "LIDPAD"],
  alternative = "two.sided", exact = FALSE, correct = FALSE
))

############################################################
#### 11. Top Gram-Stratified Bar Plots ####
############################################################

ps_trf_gram_subset <- tax_fix(ps_cqn_flfl1, min_length = 2)
ps_trf_gram_subset <- tax_transform(ps_trf_gram_subset, trans = "compositional", rank = "Genus")

ps_trf_gram_neg <- subset_taxa(ps_trf_gram_subset, tax_annot.1$GramStain %in% "Gram negative")
ps_trf_gram_pos <- subset_taxa(ps_trf_gram_subset, tax_annot.1$GramStain %in% "Gram positive")

cat("Gram negative sample sums:\n"); print(sample_sums(ps_trf_gram_neg))
cat("Gram positive sample sums:\n"); print(sample_sums(ps_trf_gram_pos))

## Gram-Negative: Top 15 (~98% coverage)
top_gneg        <- top_taxa(ps_trf_gram_neg, n = 15)
top_gneg_pruned <- prune_taxa(top_gneg, ps_trf_gram_neg)

cat("Gram-neg top15 coverage per sample:\n")
print(sample_sums(top_gneg_pruned) / sample_sums(ps_trf_gram_neg))

# Compositional version
top_gneg_pruned_comp <- transform_sample_counts(top_gneg_pruned, function(ASV) ASV / sum(ASV))

plot_bar(top_gneg_pruned_comp, x = "SampleID", fill = rank) +
  facet_wrap(~ Diet + Genotype, scales = "free_x", ncol = 5) +
  scale_fill_manual(values = c(
    microshades_palette("micro_cvd_gray"),
    microshades_palette("micro_cvd_blue"),
    microshades_palette("micro_cvd_orange"),
    microshades_palette("micro_brown")
  )) +
  ggtitle("Gram Negative (Top 15, ~98%)")

write.csv(otu_table(top_gneg_pruned), "top_gneg_otu.csv")

## Gram-Positive: Top 25 (~85% coverage)
top_gpos        <- top_taxa(ps_trf_gram_pos, n = 25)
top_gpos_pruned <- prune_taxa(top_gpos, ps_trf_gram_pos)

cat("Gram-pos top25 coverage per sample:\n")
print(sample_sums(top_gpos_pruned) / sample_sums(ps_trf_gram_pos))
cat("Mean coverage: ", mean(sample_sums(top_gpos_pruned) / sample_sums(ps_trf_gram_pos)), "\n")

top_gpos_pruned_comp <- transform_sample_counts(top_gpos_pruned, function(ASV) ASV / sum(ASV))

plot_bar(top_gpos_pruned_comp, x = "SampleID", fill = rank) +
  facet_wrap(~ Diet + Genotype, scales = "free_x", ncol = 5) +
  scale_fill_manual(values = c(
    microshades_palette("micro_blue"),
    microshades_palette("micro_brown"),
    microshades_palette("micro_cvd_gray"),
    microshades_palette("micro_cvd_green"),
    microshades_palette("micro_cvd_purple")
  )) +
  ggtitle("Gram Positive (Top 25, ~85%)")

write.csv(otu_table(top_gpos_pruned), "top_gpos_otu.csv")



############################################################
#### 12. PICRUSt2 Export & Functional Pathway Analysis ####
############################################################
# Ref: https://github.com/picrust/picrust2
# Ref: https://github.com/cafferychen777/ggpicrust2

library(phyloseq)
library(gautils)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(pheatmap)

##### 12a. Export phyloseq objects for PICRUSt2 #####

phyloseq_to_picrust2(ps_cqn_omit1, output.dir = "picrust2_input/all")


##### 12b. Load PICRUSt2 outputs & metadata #####

KO_FILE <- "./picrust2_input/picrust2_out_pipeline/bac_KO_predicted.tsv"
EC_FILE <- "./picrust2_input/picrust2_out_pipeline/bac_EC_predicted.tsv"

metadata <- as.data.frame(sample_data(ps_cqn_omit1))
metadata$Group <- paste(metadata$Genotype, metadata$Diet, sep = "_")
cat("Metadata dimensions:", dim(metadata), "\n")

# Load ASV-level KO predictions (ASVs x KOs)
ko_asv <- read_tsv(KO_FILE)
cat("PICRUSt2 output: ", nrow(ko_asv), "ASVs x", ncol(ko_asv) - 1, "KOs\n")

##### 12c. Aggregate KO predictions to sample level #####

# Phyloseq OTU table: samples x ASVs
otu <- as.matrix(otu_table(ps_cqn_omit1))
if (!taxa_are_rows(ps_cqn_omit1)) {
  otu <- t(otu)  # ensure samples x ASVs orientation
}
cat("Phyloseq OTU table: ", nrow(otu), "samples x", ncol(otu), "ASVs\n")

# Match ASVs between phyloseq and PICRUSt2
common_asvs <- intersect(colnames(otu), ko_asv$sequence)
cat("Matched ASVs:", length(common_asvs), "\n")

otu_subset <- otu[, common_asvs, drop = FALSE]
ko_subset  <- as.matrix(ko_asv[ko_asv$sequence %in% common_asvs, -1, drop = FALSE])

# Matrix multiply: (samples x ASVs) %*% (ASVs x KOs) -> samples x KOs
sample_kos <- otu_subset %*% ko_subset

# Format for ggpicrust2: KOs as rows, samples as columns
ko_sample_final <- data.frame(
  `#NAME` = paste0("ko:", colnames(ko_subset)),
  t(sample_kos),
  check.names = FALSE
)

cat("Aggregated KO table: ", nrow(ko_sample_final), "KOs x",
    ncol(ko_sample_final) - 1, "samples\n")

# Align metadata rows to PICRUSt2 sample columns
picrust_samples <- colnames(ko_sample_final)[-1]
metadata_aligned <- metadata[picrust_samples, , drop = FALSE]

##### 12d. Differential Abundance Analysis with ggpicrust2 #####

results <- ggpicrust2(
  data      = ko_sample_final,
  metadata  = metadata_aligned,
  group     = "Group",
  pathway   = "KO",
  daa_method         = "DESeq2",
  p.adjust           = "BH",
  p_values_threshold = 1
)

# Inspect top result
example_plot    <- results[[1]]$plot
example_results <- results[[1]]$results
print(example_plot)

##### 12e. Heatmap: Top Up/Down KOs (flfl Ctrl vs LIDPAD) #####

meta_df      <- as.data.frame(metadata)
flfl_samples <- rownames(meta_df)[meta_df$Genotype == "Angptl4.flfl"]
meta_flfl    <- meta_df[flfl_samples, , drop = FALSE]

daa_df <- results$daa_results_df %>%
  filter(
    group1 == "Angptl4.flfl_Ctrl",
    group2 == "Angptl4.flfl_LIDPAD",
    !is.na(p_values),
    !is.na(log2_fold_change)
  )

# Top 10 upregulated in LIDPAD (highest positive log2FC)
top10_up <- daa_df %>%
  arrange(desc(log2_fold_change), p_adjust) %>%
  slice_head(n = 10) %>%
  pull(feature)

# Top 10 downregulated in LIDPAD (lowest/most-negative log2FC)
top10_down <- daa_df %>%
  arrange(log2_fold_change, p_adjust) %>%
  slice_head(n = 10) %>%
  pull(feature)

all_features <- c(top10_up, top10_down)
cat("Top upregulated KOs:  ", paste(head(top10_up,  3), collapse = ", "), "\n")
cat("Top downregulated KOs:", paste(head(top10_down, 3), collapse = ", "), "\n")

# Subset and log-transform abundance
abund_mat      <- as.matrix(results$abundance)
abund_selected <- abund_mat[all_features, flfl_samples, drop = FALSE]
abund_log      <- log10(abund_selected + 1)

# Column annotation
annotation_col <- data.frame(
  Diet  = meta_flfl[flfl_samples, "Diet"],
  Group = meta_flfl[flfl_samples, "Group"],
  row.names = flfl_samples
)

ann_colors <- list(
  Diet = c(Ctrl = "black", LIDPAD = "red")
)

pheatmap(
  abund_log,
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row      = 9,
  fontsize_col      = 9,
  main              = "Top 10 Up + Down KOs\n(Angptl4.flfl Ctrl vs LIDPAD)"
)

##### 12f. ggpicrust2 with KEGG pathway mapping (file-based, LinDA) #####

results_kegg <- ggpicrust2(
  file      = KO_FILE,
  metadata  = metadata_aligned,
  group     = "Group",
  pathway   = "KO",
  daa_method   = "LinDA",
  ko_to_kegg   = TRUE,
  order        = "pathway_class",
  p_values_bar = TRUE,
  x_lab        = "pathway_name"
)

print(results_kegg[[1]]$plot)
print(results_kegg[[1]]$results)


############################################################
# End of script
############################################################