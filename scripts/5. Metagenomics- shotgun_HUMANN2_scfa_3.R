#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel", "ggpubr", "viridis")
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(ggpubr)
  library(viridis)
})

# ------------------ args ------------------
args <- commandArgs(trailingOnly = TRUE)
humann_out_dir <- ifelse(length(args) >= 1, args[1], ".")
out_dir        <- ifelse(length(args) >= 2, args[2], "humann_plots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------ helpers ------------------
safe_pdf_device <- function(filename, width = 10, height = 8) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename, width = width, height = height)
  } else {
    grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  }
}

safe_ggsave_pdf <- function(plot, filename, width = 10, height = 8) {
  if (capabilities("cairo")) {
    ggplot2::ggsave(filename, plot = plot, device = grDevices::cairo_pdf, width = width, height = height)
  } else {
    ggplot2::ggsave(filename, plot = plot, device = "pdf", width = width, height = height)
  }
}

normalize_pathway_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\|.*$", "", x)       # drop anything after |
  x <- sub(":.*$", "", x)         # drop anything after :
  x <- sub("\\s+.*$", "", x)      # drop anything after whitespace
  x <- trimws(x)
  toupper(x)
}

# ------------------ inputs ------------------
meta_file <- file.path(humann_out_dir, "humann_metadata.tsv")

candidate_path_files <- c(
  file.path(humann_out_dir, "humann_pathabundance_cpm_unstratified.tsv"),
  file.path(humann_out_dir, "humann_pathabundance_cpm.tsv"),
  file.path(humann_out_dir, "humann_pathabundance.tsv"),
  file.path(humann_out_dir, "pathabundance.tsv")
)
path_file <- candidate_path_files[file.exists(candidate_path_files)][1]
if (is.na(path_file) || !file.exists(path_file)) {
  stop("Could not find a pathway abundance file. Tried:\n  ", paste(candidate_path_files, collapse = "\n  "))
}
if (!file.exists(meta_file)) stop("Missing metadata: ", meta_file)

message("Using pathway file: ", path_file)
message("Using metadata file: ", meta_file)
message("Outputs -> ", normalizePath(out_dir))

# ------------------ metadata ------------------
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)
if (!all(c("sample_id", "group") %in% colnames(meta))) {
  stop("Metadata must have columns: sample_id, group")
}

inu_groups <- c("Inu_Ctrl", "Inu_LP", "Inu_LPI")
meta_inu <- meta %>%
  mutate(group = as.character(group)) %>%
  filter(group %in% inu_groups) %>%
  mutate(group = factor(group, levels = inu_groups))

if (nrow(meta_inu) < 2) stop("Too few Inu samples in metadata for groups: ", paste(inu_groups, collapse = ", "))

# ------------------ read pathway table (UNSTRATIFIED only) ------------------
pathway_df <- readr::read_tsv(path_file, show_col_types = FALSE)
path_col <- colnames(pathway_df)[1]

pathway_df_clean <- pathway_df %>%
  filter(
    !str_detect(.data[[path_col]], "UNMAPPED|UNINTEGRATED|UNGROUPED"),
    !str_detect(.data[[path_col]], "\\|")
  ) %>%
  mutate(
    pathway_full = .data[[path_col]],
    pathway_id   = normalize_pathway_id(.data[[path_col]])
  )

sample_cols <- setdiff(colnames(pathway_df_clean), c(path_col, "pathway_full", "pathway_id"))
common_samples <- intersect(sample_cols, meta_inu$sample_id)
if (length(common_samples) < 2) stop("Fewer than 2 overlapping Inu samples between pathway table and metadata.")

pathway_df_clean <- pathway_df_clean %>%
  select(pathway_full, pathway_id, all_of(common_samples))

meta_inu <- meta_inu %>% filter(sample_id %in% common_samples)
meta_inu <- meta_inu[match(common_samples, meta_inu$sample_id), ]
meta_inu$group <- factor(meta_inu$group, levels = inu_groups)

# Collapse duplicates
path_mat <- pathway_df_clean %>%
  select(pathway_id, all_of(common_samples)) %>%
  group_by(pathway_id) %>%
  summarise(across(all_of(common_samples), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("pathway_id") %>%
  as.matrix()

denom <- colSums(path_mat, na.rm = TRUE)
denom[denom == 0] <- NA_real_

# ------------------ SCFA/Butyrate pathway identification ------------------
scfa_keywords <- c(
  "butyrate", "butanoate", "propanoate", "propionate",
  "acetate", "acetyl-coa", "short chain fatty acid",
  "fermentation", "succinate", "propanediol"
)

scfa_ids <- toupper(c(
  "PWY-5022", "PWY-5676", "PWY-5677", "CENTFERM-PWY", "PWY-7013",
  "P162-PWY", "P461-PWY", "FERMENTATION-PWY"
))

id_to_full <- pathway_df_clean %>%
  distinct(pathway_id, pathway_full) %>%
  group_by(pathway_id) %>%
  summarise(pathway_full = pathway_full[1], .groups = "drop")

all_ids <- rownames(path_mat)
full_map <- id_to_full$pathway_full
names(full_map) <- id_to_full$pathway_id
path_full_for_id <- full_map[all_ids]
path_full_for_id[is.na(path_full_for_id)] <- all_ids[is.na(path_full_for_id)]

kw_regex <- paste(scfa_keywords, collapse = "|")
kw_hit <- str_detect(tolower(path_full_for_id), tolower(kw_regex))
id_hit <- all_ids %in% scfa_ids

scfa_targets <- all_ids[kw_hit | id_hit]
scfa_targets <- unique(scfa_targets)

if (length(scfa_targets) == 0) {
  stop("No SCFA-related pathways found by keyword/ID matching.")
}

message("Found ", length(scfa_targets), " SCFA-related pathways (Inu-only table).")

# ------------------ normalization ------------------
scfa_mat <- path_mat[scfa_targets, , drop = FALSE]
scfa_rel <- sweep(scfa_mat, 2, denom, FUN = "/") * 100
scfa_log <- log10(scfa_mat + 1)

# ------------------ statistics (3-group) ------------------
group_vec <- meta_inu$group
names(group_vec) <- meta_inu$sample_id

kw_p <- apply(scfa_rel, 1, function(x) {
  x <- as.numeric(x)
  suppressWarnings(stats::kruskal.test(x ~ group_vec)$p.value)
})
kw_fdr <- p.adjust(kw_p, method = "BH")

pairwise_for_pathway <- function(x) {
  x <- as.numeric(x)
  pw <- suppressWarnings(stats::pairwise.wilcox.test(x, group_vec, p.adjust.method = "BH"))
  out <- c(
    p_LP_vs_Ctrl  = pw$p.value["Inu_LP",  "Inu_Ctrl"],
    p_LPI_vs_LP   = pw$p.value["Inu_LPI", "Inu_LP"],
    p_LPI_vs_Ctrl = pw$p.value["Inu_LPI", "Inu_Ctrl"]
  )
  as.numeric(out)
}

pw_mat <- t(apply(scfa_rel, 1, pairwise_for_pathway))
colnames(pw_mat) <- c("p_LP_vs_Ctrl_BH", "p_LPI_vs_LP_BH", "p_LPI_vs_Ctrl_BH")

means_by_group <- function(mat) {
  tibble(
    pathway_id = rownames(mat),
    Mean_Inu_Ctrl = rowMeans(mat[, meta_inu$sample_id[meta_inu$group == "Inu_Ctrl"], drop = FALSE], na.rm = TRUE),
    Mean_Inu_LP   = rowMeans(mat[, meta_inu$sample_id[meta_inu$group == "Inu_LP"],   drop = FALSE], na.rm = TRUE),
    Mean_Inu_LPI  = rowMeans(mat[, meta_inu$sample_id[meta_inu$group == "Inu_LPI"],  drop = FALSE], na.rm = TRUE)
  )
}

means_rel <- means_by_group(scfa_rel)

stat_df <- means_rel %>%
  mutate(
    KW_p = as.numeric(kw_p[pathway_id]),
    KW_FDR = as.numeric(kw_fdr[pathway_id])
  ) %>%
  bind_cols(as.data.frame(pw_mat)) %>%
  mutate(
    Fold_LP_vs_Ctrl  = (Mean_Inu_LP  + 0.01) / (Mean_Inu_Ctrl + 0.01),
    Fold_LPI_vs_Ctrl = (Mean_Inu_LPI + 0.01) / (Mean_Inu_Ctrl + 0.01),
    Fold_LPI_vs_LP   = (Mean_Inu_LPI + 0.01) / (Mean_Inu_LP   + 0.01),
    pathway_full = path_full_for_id[pathway_id]
  ) %>%
  relocate(pathway_id, pathway_full) %>%
  arrange(KW_p)

stat_csv <- file.path(out_dir, "SCFA_pathway_statistics_InuOnly.csv")
write.csv(stat_df, stat_csv, row.names = FALSE)
message("Saved: ", stat_csv)

# ------------------ ORIGINAL PLOTS ------------------
annotation_col <- data.frame(
  Group = meta_inu$group,
  row.names = meta_inu$sample_id
)
annotation_colors <- list(
  Group = c(Inu_Ctrl = "#56B4E9", Inu_LP = "#E69F00", Inu_LPI = "#009E73")
)

# HEATMAP - with robust data cleaning
pathway_rowsums <- rowSums(scfa_mat, na.rm = TRUE)
pathway_variance <- apply(scfa_mat, 1, var, na.rm = TRUE)
valid_pathways <- (pathway_rowsums > 0) & (!is.na(pathway_variance)) & (pathway_variance > 0)

scfa_log_clean <- scfa_log[valid_pathways, , drop = FALSE]

scfa_log_clean[is.na(scfa_log_clean) | is.nan(scfa_log_clean) | is.infinite(scfa_log_clean)] <- 0

if (nrow(scfa_log_clean) >= 2 && ncol(scfa_log_clean) >= 2) {
  heat_pdf <- file.path(out_dir, "SCFA_pathways_heatmap_InuOnly.pdf")
  safe_pdf_device(heat_pdf, width = 10, height = 8)
  pheatmap(
    scfa_log_clean,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    main = "SCFA/Butyrate Pathways (log10 CPM+1; row-scaled) | Inu groups only",
    fontsize = 10,
    fontsize_row = 8
  )
  dev.off()
  message("Saved: ", heat_pdf)
} else {
  message("Skipping heatmap: insufficient valid pathways for clustering")
}


top_n <- min(6, nrow(stat_df))
top_pathways <- stat_df$pathway_id[seq_len(top_n)]

plot_df <- scfa_rel[top_pathways, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("pathway_id") %>%
  pivot_longer(-pathway_id, names_to = "sample_id", values_to = "abundance_rel") %>%
  left_join(meta_inu, by = "sample_id") %>%
  mutate(
    pathway_id = factor(pathway_id, levels = top_pathways),
    group = factor(group, levels = inu_groups)
  )

box_pdf <- file.path(out_dir, "SCFA_pathways_boxplots_top6_InuOnly.pdf")
p_box <- ggplot(plot_df, aes(x = group, y = abundance_rel, fill = group)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.18, size = 1.6, alpha = 0.55) +
  facet_wrap(~ pathway_id, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = annotation_colors$Group) +
  labs(
    title = "Top SCFA/Butyrate pathways (by Kruskal-Wallis p) | Inu groups",
    x = "Group", y = "Relative abundance (%)", fill = "Group"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

safe_ggsave_pdf(p_box, box_pdf, width = 12, height = 10)
message("Saved: ", box_pdf)

# VOLCANO PLOTS
make_volcano <- function(df_stats, comparison = c("LP_vs_Ctrl", "LPI_vs_Ctrl"), out_pdf) {
  comparison <- match.arg(comparison)
  
  if (comparison == "LP_vs_Ctrl") {
    log2fc <- log2(df_stats$Fold_LP_vs_Ctrl)
    pcol   <- "p_LP_vs_Ctrl_BH"
    title  <- "Volcano (Inu_LP vs Inu_Ctrl)"
  } else {
    log2fc <- log2(df_stats$Fold_LPI_vs_Ctrl)
    pcol   <- "p_LPI_vs_Ctrl_BH"
    title  <- "Volcano (Inu_LPI vs Inu_Ctrl)"
  }
  
  v <- df_stats %>%
    mutate(
      log2FC = log2fc,
      p_adj = .data[[pcol]],
      neglog10p = -log10(p_adj + 1e-300),
      Category = case_when(
        p_adj < 0.05 & log2FC >  0.5 ~ "Enriched (treatment)",
        p_adj < 0.05 & log2FC < -0.5 ~ "Enriched (Ctrl)",
        TRUE ~ "Not significant"
      ),
      label = ifelse(p_adj < 0.05 & abs(log2FC) > 0.5, pathway_id, "")
    )
  
  p <- ggplot(v, aes(x = log2FC, y = neglog10p, color = Category)) +
    geom_point(size = 2.6, alpha = 0.75) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.7) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 2.6,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.25,
      force = 4,
      min.segment.length = 0,
      segment.alpha = 0.35,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c(
      "Enriched (treatment)" = "#E69F00",
      "Enriched (Ctrl)"      = "#56B4E9",
      "Not significant"      = "grey60"
    )) +
    labs(
      title = title,
      x = "log2(Fold change)",
      y = "-log10(adj p)",
      color = "Category"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  safe_ggsave_pdf(p, out_pdf, width = 10, height = 8)
  message("Saved: ", out_pdf)
}

vol1 <- file.path(out_dir, "SCFA_pathways_volcano_Inu_LP_vs_Inu_Ctrl.pdf")
vol2 <- file.path(out_dir, "SCFA_pathways_volcano_Inu_LPI_vs_Inu_Ctrl.pdf")
make_volcano(stat_df, "LP_vs_Ctrl",  vol1)
make_volcano(stat_df, "LPI_vs_Ctrl", vol2)

# TOTAL SCFA BARPLOT
total_scfa <- tibble(
  sample_id = colnames(scfa_rel),
  Total_SCFA_rel = colSums(scfa_rel, na.rm = TRUE)
) %>%
  left_join(meta_inu, by = "sample_id") %>%
  mutate(group = factor(group, levels = inu_groups))

bar_pdf <- file.path(out_dir, "Total_SCFA_abundance_barplot_InuOnly.pdf")
p_bar <- ggplot(total_scfa, aes(x = reorder(sample_id, Total_SCFA_rel), y = Total_SCFA_rel, fill = group)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = annotation_colors$Group) +
  labs(
    title = "Total SCFA/Butyrate pathway abundance per sample (sum of selected pathways)",
    x = "Sample", y = "Total relative abundance (%)", fill = "Group"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")

safe_ggsave_pdf(p_bar, bar_pdf, width = 12, height = 6)
message("Saved: ", bar_pdf)


message("\n=== Generating plots ===")

# ------------------ PLOT 1: Group comparison of total SCFA capacity ------------------

total_scfa_stats <- total_scfa %>%
  group_by(group) %>%
  summarise(
    mean = mean(Total_SCFA_rel, na.rm = TRUE),
    sd = sd(Total_SCFA_rel, na.rm = TRUE),
    se = sd / sqrt(n()),
    .groups = "drop"
  )

# Statistical test for total SCFA
kw_total <- kruskal.test(Total_SCFA_rel ~ group, data = total_scfa)
pw_total <- pairwise.wilcox.test(total_scfa$Total_SCFA_rel, total_scfa$group, p.adjust.method = "BH")

p_total_comp <- ggplot(total_scfa, aes(x = group, y = Total_SCFA_rel, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white", color = "black") +
  scale_fill_manual(values = annotation_colors$Group) +
  stat_compare_means(method = "kruskal.test", label.y.npc = 0.95, size = 4.5) +
  stat_compare_means(
    comparisons = list(c("Inu_Ctrl", "Inu_LP"), c("Inu_LP", "Inu_LPI"), c("Inu_Ctrl", "Inu_LPI")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.08
  ) +
  labs(
    title = "Total SCFA Pathway Capacity Across Groups",
    subtitle = "Diamond = group mean; box = median and IQR",
    x = "Treatment Group",
    y = "Total SCFA Pathway Abundance (%)",
    fill = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11)
  )

total_comp_pdf <- file.path(out_dir, "Total_SCFA_group_comparison.pdf")
safe_ggsave_pdf(p_total_comp, total_comp_pdf, width = 8, height = 7)
message("Saved: ", total_comp_pdf)

# ------------------ PLOT 2: Mean abundance profiles across groups ------------------

mean_profiles <- stat_df %>%
  select(pathway_id, Mean_Inu_Ctrl, Mean_Inu_LP, Mean_Inu_LPI) %>%
  pivot_longer(-pathway_id, names_to = "group", values_to = "mean_abundance") %>%
  mutate(
    group = factor(group, 
                   levels = c("Mean_Inu_Ctrl", "Mean_Inu_LP", "Mean_Inu_LPI"),
                   labels = c("Inu_Ctrl", "Inu_LP", "Inu_LPI"))
  )

# Identify top pathways by max mean abundance
top_by_abundance <- stat_df %>%
  mutate(max_mean = pmax(Mean_Inu_Ctrl, Mean_Inu_LP, Mean_Inu_LPI)) %>%
  arrange(desc(max_mean)) %>%
  head(10) %>%
  pull(pathway_id)

profile_data <- mean_profiles %>%
  filter(pathway_id %in% top_by_abundance) %>%
  mutate(pathway_id = factor(pathway_id, levels = top_by_abundance))

p_profiles <- ggplot(profile_data, aes(x = pathway_id, y = mean_abundance, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.85) +
  scale_fill_manual(values = annotation_colors$Group) +
  labs(
    title = "Top 10 SCFA Pathways: Mean Abundance Profiles",
    subtitle = "Which pathways dominate in each group?",
    x = "Pathway ID",
    y = "Mean Relative Abundance (%)",
    fill = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

profiles_pdf <- file.path(out_dir, "SCFA_pathways_mean_profiles_top10.pdf")
safe_ggsave_pdf(p_profiles, profiles_pdf, width = 12, height = 7)
message("Saved: ", profiles_pdf)

# ------------------ PLOT 3: Pathway enrichment summary heatmap ------------------

fc_data <- stat_df %>%
  select(pathway_id, Fold_LP_vs_Ctrl, Fold_LPI_vs_Ctrl, Fold_LPI_vs_LP,
         p_LP_vs_Ctrl_BH, p_LPI_vs_Ctrl_BH, p_LPI_vs_LP_BH) %>%
  mutate(
    log2FC_LP_vs_Ctrl = log2(Fold_LP_vs_Ctrl),
    log2FC_LPI_vs_Ctrl = log2(Fold_LPI_vs_Ctrl),
    log2FC_LPI_vs_LP = log2(Fold_LPI_vs_LP)
  )

# Create matrix for heatmap
fc_mat <- fc_data %>%
  select(pathway_id, log2FC_LP_vs_Ctrl, log2FC_LPI_vs_Ctrl, log2FC_LPI_vs_LP) %>%
  column_to_rownames("pathway_id") %>%
  as.matrix()

colnames(fc_mat) <- c("LP vs Ctrl", "LPI vs Ctrl", "LPI vs LP")

# Clean infinite and NA values
fc_mat[is.infinite(fc_mat)] <- NA
fc_mat[is.na(fc_mat)] <- 0  # Replace NA with 0 (no change)

# Create significance annotation
sig_mat <- fc_data %>%
  select(pathway_id, p_LP_vs_Ctrl_BH, p_LPI_vs_Ctrl_BH, p_LPI_vs_LP_BH) %>%
  column_to_rownames("pathway_id") %>%
  as.matrix()

sig_mat_display <- sig_mat
sig_mat_display[sig_mat < 0.001] <- "***"
sig_mat_display[sig_mat >= 0.001 & sig_mat < 0.01] <- "**"
sig_mat_display[sig_mat >= 0.01 & sig_mat < 0.05] <- "*"
sig_mat_display[sig_mat >= 0.05 | is.na(sig_mat)] <- ""

if (nrow(fc_mat) >= 2) {
  fc_heat_pdf <- file.path(out_dir, "SCFA_pathways_fold_change_heatmap.pdf")
  safe_pdf_device(fc_heat_pdf, width = 8, height = max(6, nrow(fc_mat) * 0.3))
  pheatmap(
    fc_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-2, 2, length.out = 101),
    main = "SCFA Pathway Fold Changes (log2)\nAcross All Comparisons",
    fontsize = 9,
    fontsize_row = 7,
    display_numbers = sig_mat_display,
    number_color = "black",
    fontsize_number = 8,
    cellwidth = 40,
    cellheight = 12
  )
  dev.off()
  message("Saved: ", fc_heat_pdf)
} else {
  message("Skipping fold change heatmap: insufficient pathways")
}

# ------------------ PLOT 4: Progressive treatment effect plot ------------------

# Calculate effect sizes
effect_sizes <- tibble(
  Comparison = c("Baseline", "Inu_LP (+probiotic)", "Inu_LPI (+synbiotic)"),
  Mean_Total_SCFA = c(
    mean(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_Ctrl"]),
    mean(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_LP"]),
    mean(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_LPI"])
  ),
  SE = c(
    sd(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_Ctrl"]) / sqrt(sum(total_scfa$group == "Inu_Ctrl")),
    sd(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_LP"]) / sqrt(sum(total_scfa$group == "Inu_LP")),
    sd(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_LPI"]) / sqrt(sum(total_scfa$group == "Inu_LPI"))
  )
) %>%
  mutate(
    Comparison = factor(Comparison, levels = c("Baseline", "Inu_LP (+probiotic)", "Inu_LPI (+synbiotic)")),
    Percent_Change = (Mean_Total_SCFA / Mean_Total_SCFA[1] - 1) * 100
  )

p_progressive <- ggplot(effect_sizes, aes(x = Comparison, y = Mean_Total_SCFA, group = 1)) +
  geom_line(linewidth = 1.2, color = "#E69F00") +
  geom_point(size = 5, color = "#E69F00") +
  geom_errorbar(aes(ymin = Mean_Total_SCFA - SE, ymax = Mean_Total_SCFA + SE), 
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percent_Change)), 
            vjust = -1.5, size = 4, fontface = "bold") +
  labs(
    title = "Progressive Effect of Probiotic Addition on SCFA Capacity",
    subtitle = "How does SCFA production change from baseline through interventions?",
    x = "Treatment Progression",
    y = "Mean Total SCFA Pathway Abundance (%) +/- SE",
    caption = "Numbers show percent change from baseline"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11)
  ) +
  ylim(NA, max(effect_sizes$Mean_Total_SCFA + effect_sizes$SE) * 1.15)

progressive_pdf <- file.path(out_dir, "SCFA_progressive_treatment_effect.pdf")
safe_ggsave_pdf(p_progressive, progressive_pdf, width = 9, height = 7)
message("Saved: ", progressive_pdf)

# ------------------ PLOT 5: Individual pathway response patterns ------------------

pathway_patterns <- stat_df %>%
  mutate(
    Pattern = case_when(
      # Both treatments increase
      Fold_LP_vs_Ctrl > 1.2 & Fold_LPI_vs_Ctrl > 1.2 ~ "Increased by both LP and LPI",
      # Only LP increases
      Fold_LP_vs_Ctrl > 1.2 & Fold_LPI_vs_Ctrl <= 1.2 ~ "Increased by LP only",
      # Only LPI increases
      Fold_LP_vs_Ctrl <= 1.2 & Fold_LPI_vs_Ctrl > 1.2 ~ "Increased by LPI only",
      # LPI additive effect (higher than LP)
      Fold_LP_vs_Ctrl > 1.2 & Fold_LPI_vs_LP > 1.2 ~ "LPI shows additive effect",
      # Both decrease
      Fold_LP_vs_Ctrl < 0.8 & Fold_LPI_vs_Ctrl < 0.8 ~ "Decreased by both",
      # No clear pattern
      TRUE ~ "No consistent change"
    )
  )

pattern_summary <- pathway_patterns %>%
  count(Pattern) %>%
  arrange(desc(n))

p_patterns <- ggplot(pattern_summary, aes(x = reorder(Pattern, n), y = n, fill = Pattern)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = n), hjust = -0.2, size = 4.5, fontface = "bold") +
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "SCFA Pathway Response Patterns",
    subtitle = "How do individual pathways respond to treatments?",
    x = "Response Pattern",
    y = "Number of Pathways",
    fill = "Pattern"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  expand_limits(y = max(pattern_summary$n) * 1.15)

patterns_pdf <- file.path(out_dir, "SCFA_pathway_response_patterns.pdf")
safe_ggsave_pdf(p_patterns, patterns_pdf, width = 10, height = 6)
message("Saved: ", patterns_pdf)

# Save pattern classifications
pattern_csv <- file.path(out_dir, "SCFA_pathway_patterns.csv")
write.csv(pathway_patterns, pattern_csv, row.names = FALSE)
message("Saved: ", pattern_csv)

# ------------------ PLOT 6: Specific SCFA type breakdown  ------------------

scfa_type_classification <- stat_df %>%
  mutate(
    SCFA_Type = case_when(
      str_detect(tolower(pathway_full), "butyrate|butanoate") ~ "Butyrate",
      str_detect(tolower(pathway_full), "propionate|propanoate") ~ "Propionate",
      str_detect(tolower(pathway_full), "acetate|acetyl") ~ "Acetate",
      TRUE ~ "Mixed/Other"
    )
  )

# Calculate total abundance by SCFA type for each sample
scfa_type_abundance <- scfa_rel %>%
  as.data.frame() %>%
  rownames_to_column("pathway_id") %>%
  left_join(scfa_type_classification %>% select(pathway_id, SCFA_Type), by = "pathway_id") %>%
  pivot_longer(-c(pathway_id, SCFA_Type), names_to = "sample_id", values_to = "abundance") %>%
  group_by(SCFA_Type, sample_id) %>%
  summarise(Total_Abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  left_join(meta_inu, by = "sample_id")

# Summary stats by type and group
scfa_type_summary <- scfa_type_abundance %>%
  group_by(SCFA_Type, group) %>%
  summarise(
    Mean = mean(Total_Abundance, na.rm = TRUE),
    SE = sd(Total_Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_scfa_types <- ggplot(scfa_type_summary, aes(x = group, y = Mean, fill = SCFA_Type)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "SCFA Production by Type Across Groups",
    subtitle = "Breakdown of butyrate, propionate, and acetate pathway capacity",
    x = "Treatment Group",
    y = "Mean Total Pathway Abundance (%) +/- SE",
    fill = "SCFA Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11)
  )

scfa_types_pdf <- file.path(out_dir, "SCFA_production_by_type.pdf")
safe_ggsave_pdf(p_scfa_types, scfa_types_pdf, width = 10, height = 7)
message("Saved: ", scfa_types_pdf)

# ------------------ PLOT 7: Sample-level heatmap with total SCFA annotation ------------------

# Prepare annotation with total SCFA
annotation_col_enhanced <- data.frame(
  Group = meta_inu$group,
  Total_SCFA = total_scfa$Total_SCFA_rel[match(meta_inu$sample_id, total_scfa$sample_id)],
  row.names = meta_inu$sample_id
)

annotation_colors_enhanced <- list(
  Group = c(Inu_Ctrl = "#56B4E9", Inu_LP = "#E69F00", Inu_LPI = "#009E73"),
  Total_SCFA = viridis(100)
)

# Filter to top pathways by variance for cleaner visualization
pathway_variance <- apply(scfa_log_clean, 1, var, na.rm = TRUE)
pathway_variance[is.na(pathway_variance)] <- 0
top_var_pathways <- names(sort(pathway_variance, decreasing = TRUE))[1:min(20, length(pathway_variance))]
top_var_pathways <- top_var_pathways[!is.na(top_var_pathways)]

if (length(top_var_pathways) >= 2) {
  scfa_log_top <- scfa_log_clean[top_var_pathways, , drop = FALSE]
  
  heat_enhanced_pdf <- file.path(out_dir, "SCFA_pathways_heatmap_enhanced.pdf")
  safe_pdf_device(heat_enhanced_pdf, width = 11, height = 9)
  pheatmap(
    scfa_log_top,
    annotation_col = annotation_col_enhanced,
    annotation_colors = annotation_colors_enhanced,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    main = "Top Variable SCFA Pathways\n(annotated by group and total SCFA capacity)",
    fontsize = 10,
    fontsize_row = 8
  )
  dev.off()
  message("Saved: ", heat_enhanced_pdf)
} else {
  message("Skipping enhanced heatmap: insufficient variable pathways")
}

# ------------------ PLOT 8: Correlation between pathways ------------------

if (nrow(scfa_rel) >= 3 && ncol(scfa_rel) >= 3) {
  # Use cleaned data for correlation
  scfa_rel_clean <- scfa_rel[valid_pathways, , drop = FALSE]
  
  # Only calculate correlation if we have enough data
  if (nrow(scfa_rel_clean) >= 3) {
    cor_mat <- cor(t(scfa_rel_clean), use = "pairwise.complete.obs", method = "spearman")
    
    # Remove rows/cols that are all NA
    valid_cor <- !is.na(rowSums(cor_mat))
    cor_mat <- cor_mat[valid_cor, valid_cor, drop = FALSE]
    
    # Filter for visualization if too many pathways
    if (nrow(cor_mat) > 30) {
      top_pathways_cor <- rownames(cor_mat)[1:min(30, nrow(cor_mat))]
      cor_mat <- cor_mat[top_pathways_cor, top_pathways_cor, drop = FALSE]
    }
    
    if (nrow(cor_mat) >= 2) {
      cor_pdf <- file.path(out_dir, "SCFA_pathway_correlations.pdf")
      safe_pdf_device(cor_pdf, width = 10, height = 9)
      pheatmap(
        cor_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        breaks = seq(-1, 1, length.out = 101),
        main = "SCFA Pathway Co-occurrence Patterns\n(Spearman correlation)",
        fontsize_row = 7,
        fontsize_col = 7
      )
      dev.off()
      message("Saved: ", cor_pdf)
    } else {
      message("Skipping correlation heatmap: insufficient valid pathways after cleaning")
    }
  } else {
    message("Skipping correlation heatmap: insufficient pathways")
  }
} else {
  message("Skipping correlation heatmap: insufficient data dimensions")
}

# ------------------ PLOT 9: Effect size comparison plot ------------------

effect_comparison <- stat_df %>%
  select(pathway_id, Fold_LP_vs_Ctrl, Fold_LPI_vs_Ctrl, Fold_LPI_vs_LP,
         p_LP_vs_Ctrl_BH, p_LPI_vs_Ctrl_BH, p_LPI_vs_LP_BH) %>%
  pivot_longer(
    cols = c(Fold_LP_vs_Ctrl, Fold_LPI_vs_Ctrl, Fold_LPI_vs_LP),
    names_to = "Comparison",
    values_to = "Fold_Change"
  ) %>%
  mutate(
    log2FC = log2(Fold_Change),
    Comparison = factor(
      Comparison,
      levels = c("Fold_LP_vs_Ctrl", "Fold_LPI_vs_Ctrl", "Fold_LPI_vs_LP"),
      labels = c("LP vs Ctrl", "LPI vs Ctrl", "LPI vs LP")
    )
  )

p_effect_comparison <- ggplot(effect_comparison, aes(x = Comparison, y = log2FC, fill = Comparison)) +
  geom_jitter(alpha = 0.5, width = 0.35) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) +
  labs(
    title = "Distribution of Effect Sizes Across Comparisons",
    subtitle = "Which treatment comparison shows the strongest SCFA pathway changes?",
    x = "Comparison",
    y = "log2(Fold Change)",
    fill = "Comparison"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

effect_comp_pdf <- file.path(out_dir, "SCFA_effect_size_comparison.pdf")
safe_ggsave_pdf(p_effect_comparison, effect_comp_pdf, width = 9, height = 7)
message("Saved: ", effect_comp_pdf)

# ------------------ PLOT 10: Summary dashboard figure ------------------

# Prepare data for summary
n_sig_LP <- sum(stat_df$p_LP_vs_Ctrl_BH < 0.05, na.rm = TRUE)
n_sig_LPI <- sum(stat_df$p_LPI_vs_Ctrl_BH < 0.05, na.rm = TRUE)
n_sig_LPI_LP <- sum(stat_df$p_LPI_vs_LP_BH < 0.05, na.rm = TRUE)

summary_stats_df <- tibble(
  Metric = c(
    "Total pathways analyzed",
    "Sig. different (LP vs Ctrl)",
    "Sig. different (LPI vs Ctrl)",
    "Sig. different (LPI vs LP)",
    "Mean fold change (LP vs Ctrl)",
    "Mean fold change (LPI vs Ctrl)"
  ),
  Value = c(
    nrow(stat_df),
    n_sig_LP,
    n_sig_LPI,
    n_sig_LPI_LP,
    round(mean(stat_df$Fold_LP_vs_Ctrl, na.rm = TRUE), 2),
    round(mean(stat_df$Fold_LPI_vs_Ctrl, na.rm = TRUE), 2)
  )
)

# Create multi-panel summary
p1_sum <- ggplot(total_scfa, aes(x = group, y = Total_SCFA_rel, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = annotation_colors$Group) +
  labs(title = "A. Total SCFA Capacity", y = "Total abundance (%)") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank())

p2_sum <- ggplot(pattern_summary, aes(x = reorder(Pattern, n), y = n, fill = Pattern)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_d() +
  labs(title = "B. Response Patterns", x = NULL, y = "Count") +
  theme_bw() +
  theme(legend.position = "none")

p3_sum <- ggplot(effect_sizes, aes(x = Comparison, y = Mean_Total_SCFA)) +
  geom_line(aes(group = 1), linewidth = 1.2, color = "#E69F00") +
  geom_point(size = 4, color = "#E69F00") +
  labs(title = "C. Progressive Effect", y = "Mean SCFA (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 8))

p4_sum <- ggplot(summary_stats_df, aes(x = reorder(Metric, Value), y = Value)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  geom_text(aes(label = Value), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(title = "D. Summary Statistics", x = NULL, y = "Value") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))

combined_summary <- ggpubr::ggarrange(
  p1_sum, p2_sum, p3_sum, p4_sum,
  ncol = 2, nrow = 2
)

summary_annotated <- ggpubr::annotate_figure(
  combined_summary,
  top = text_grob("SCFA Pathway Analysis Summary Dashboard", 
                  face = "bold", size = 16)
)

summary_pdf <- file.path(out_dir, "SCFA_analysis_summary_dashboard.pdf")
safe_ggsave_pdf(summary_annotated, summary_pdf, width = 14, height = 10)
message("Saved: ", summary_pdf)

# ------------------ Final summary statistics ------------------
message("\n=== SUMMARY (Inu-only) ===")
message("n pathways (SCFA set): ", length(scfa_targets))
message("n samples: ", nrow(meta_inu))
message("Group counts:\n", paste(capture.output(print(table(meta_inu$group))), collapse = "\n"))
message("\n--- Total SCFA Capacity ---")
message("Inu_Ctrl mean: ", round(mean(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_Ctrl"]), 2))
message("Inu_LP mean: ", round(mean(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_LP"]), 2))
message("Inu_LPI mean: ", round(mean(total_scfa$Total_SCFA_rel[total_scfa$group == "Inu_LPI"]), 2))
message("Kruskal-Wallis p-value: ", format(kw_total$p.value, scientific = TRUE, digits = 3))
message("\n--- Pairwise Comparisons (Total SCFA) ---")
message("LP vs Ctrl p-value: ", format(pw_total$p.value["Inu_LP", "Inu_Ctrl"], scientific = TRUE, digits = 3))
message("LPI vs Ctrl p-value: ", format(pw_total$p.value["Inu_LPI", "Inu_Ctrl"], scientific = TRUE, digits = 3))
message("LPI vs LP p-value: ", format(pw_total$p.value["Inu_LPI", "Inu_LP"], scientific = TRUE, digits = 3))
message("\n--- Significantly Changed Pathways ---")
message("LP vs Ctrl: ", n_sig_LP, " pathways (FDR < 0.05)")
message("LPI vs Ctrl: ", n_sig_LPI, " pathways (FDR < 0.05)")
message("LPI vs LP: ", n_sig_LPI_LP, " pathways (FDR < 0.05)")

# ------------------ save workspace ------------------
save.image(file = file.path(out_dir, "4c_humann_viz-ws.RData"))
message("Saved workspace: ", file.path(out_dir, "4c_humann_viz-ws.RData"))

message("\n=== ANALYSIS COMPLETE ===")
message("Generated ", list.files(out_dir, pattern = "\\.pdf$") %>% length(), " PDF plots")
message("All outputs written to: ", normalizePath(out_dir))
message("\nKEY BIOLOGICAL INSIGHTS:")
message("1. Check 'Total_SCFA_group_comparison.pdf' to see overall SCFA capacity differences")
message("2. Check 'SCFA_progressive_treatment_effect.pdf' to see cumulative treatment effects")
message("3. Check 'SCFA_pathway_response_patterns.pdf' to understand which pathways respond to which treatments")
message("4. Check 'SCFA_production_by_type.pdf' to see butyrate vs propionate vs acetate changes")
message("5. Check 'SCFA_analysis_summary_dashboard.pdf' for a one-page overview")