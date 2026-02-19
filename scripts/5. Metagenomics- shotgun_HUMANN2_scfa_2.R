#!/usr/bin/env Rscript

## ---------- encoding/locale safety ----------
options(stringsAsFactors = FALSE)
try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)

suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "data.table", "ggrepel", "ggalluvial")
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
  library(tidyverse)
  library(data.table)
  library(ggrepel)
  library(ggalluvial)
})

## ------------------ args ------------------
args <- commandArgs(trailingOnly = TRUE)
humann_out_dir <- ifelse(length(args) >= 1, args[1], ".")
out_dir        <- ifelse(length(args) >= 2, args[2], "humann_plots")

ws_file <- file.path(out_dir, "4a_humann_viz-ws.RData")
if (!file.exists(ws_file)) stop("Workspace not found: ", ws_file)

meta_file       <- file.path(humann_out_dir, "humann_metadata.tsv")
path_file       <- file.path(humann_out_dir, "humann_pathabundance_cpm.tsv")
stratified_file <- file.path(humann_out_dir, "humann_pathabundance_cpm_stratified.tsv")

if (!file.exists(meta_file))       stop("Missing metadata: ", meta_file)
if (!file.exists(path_file))       stop("Missing unstratified pathways: ", path_file)
if (!file.exists(stratified_file)) stop("Missing stratified pathways: ", stratified_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ------------------ helpers ------------------
safe_pdf_device <- function(filename, width = 10, height = 8) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
}

extract_genus <- function(taxon) {
  ifelse(grepl("g__", taxon),
         sub(".*g__([^\\.\\;\\|]+).*", "\\1", taxon),
         taxon)
}

is_unclassified <- function(x) {
  grepl("unclassified", x, ignore.case = TRUE)
}

is_helicobacter <- function(x) {
  grepl("^helicobacter$", trimws(tolower(as.character(x))))
}

## Normalize pathway IDs to MetaCyc-style IDs
normalize_pathway_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\|.*$", "", x)    # drop after |
  x <- sub(":.*$", "", x)      # drop after :
  x <- sub("\\s+.*$", "", x)   # drop after whitespace
  x <- trimws(x)
  toupper(x)
}

## Bootstrap mean-difference CI
boot_diff_ci <- function(x1, x2, B = 2000, seed = 123) {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 2 || length(x2) < 2) {
    return(c(diff = mean(x2) - mean(x1), lo = NA_real_, hi = NA_real_))
  }
  set.seed(seed)
  diffs <- replicate(B, mean(sample(x2, replace = TRUE)) - mean(sample(x1, replace = TRUE)))
  c(diff = mean(x2) - mean(x1),
    lo = unname(stats::quantile(diffs, 0.025, na.rm = TRUE)),
    hi = unname(stats::quantile(diffs, 0.975, na.rm = TRUE)))
}

distinct_cols <- function(n) {
  if ("hcl.colors" %in% ls("package:grDevices")) {
    grDevices::hcl.colors(n, palette = "Dark 3")
  } else {
    grDevices::rainbow(n)
  }
}

detect_scfa_pathways <- function(pathway_full_vec, pathway_id_vec) {
  txt <- paste0(tolower(pathway_id_vec), " ", tolower(pathway_full_vec))

  kw <- c(
    "butyr", "butano", "butan", "croton",
    "propion", "propan", "acrylate",
    "acetat", "acetyl", "ethano",
    "valerat", "isobut", "isoval",
    "scfa", "short-chain", "short chain",
    "ferment", "lactate", "succinate", "formate",
    "pyruvate", "acetogenesis"
  )

  keep <- rep(FALSE, length(txt))
  for (k in kw) keep <- keep | grepl(k, txt, fixed = TRUE)

  unique(toupper(pathway_id_vec[keep]))
}

## ------------------ load workspace ------------------
message("Loading workspace: ", ws_file)
load(ws_file)

## ------------------ metadata ------------------
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)
if (!all(c("sample_id", "group") %in% colnames(meta))) stop("Metadata must have columns: sample_id, group")

group_order <- c("2w_Ctrl","2w_LP","Inu_Ctrl","Inu_LP","Inu_LPI","Rev_Ctrl","Rev_LP","Rev_R")
meta$group <- factor(meta$group, levels = group_order)

inu_groups <- c("Inu_Ctrl", "Inu_LP", "Inu_LPI")

## ------------------ target pathways (IDs) ------------------
target_pathways <- c(
  "PWY-6695","PWY-6151","PWY-5005","P41-PWY","GLYCOLYSIS","PWY-5484","FERMENTATION-PWY",
  "NAGPWY","PWY-5981","GALACT-GLUCUROCAT-PWY","COBALSYN-PWY","PWY-6143","HEMESYN2-PWY",
  "PWY-5188","PWY-5189","PWY-5100","OANTIGEN-PWY","ARGSYN-PWY","ARGSYNBSUB-PWY",
  "GLUTORN-PWY","PWY-5667","ASSIMI-241","PWY-5941"
)
target_pathways <- toupper(target_pathways)

## ------------------ read unstratified table (robust matching) ------------------
path_df <- readr::read_tsv(path_file, show_col_types = FALSE)
path_col <- colnames(path_df)[1]

path_df_clean <- path_df %>%
  filter(!str_detect(.data[[path_col]], "UNMAPPED|UNINTEGRATED|UNGROUPED"),
         !str_detect(.data[[path_col]], "\\|")) %>%
  mutate(
    pathway_full = .data[[path_col]],
    pathway_id   = normalize_pathway_id(.data[[path_col]])
  )

sample_cols <- setdiff(colnames(path_df_clean), c(path_col, "pathway_full", "pathway_id"))
common_samples <- intersect(sample_cols, meta$sample_id)
if (length(common_samples) < 2) stop("Fewer than 2 overlapping samples between pathways and metadata.")

path_df_clean <- path_df_clean %>%
  select(pathway_full, pathway_id, all_of(common_samples))

meta2 <- meta %>% filter(sample_id %in% common_samples)
meta2 <- meta2[match(common_samples, meta2$sample_id), ]
meta2$group <- factor(meta2$group, levels = group_order)

## Collapse duplicates at pathway_id (sum)
path_mat_all <- path_df_clean %>%
  select(pathway_id, all_of(common_samples)) %>%
  group_by(pathway_id) %>%
  summarise(across(all_of(common_samples), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("pathway_id") %>%
  as.matrix()

meta_inu <- meta2 %>% filter(group %in% inu_groups)
if (nrow(meta_inu) < 2) stop("Too few Inu samples after overlap filtering.")

## present targets
if (is.null(target_pathways)) {
  present_targets <- rownames(path_mat_all)
} else {
  present_targets <- intersect(target_pathways, rownames(path_mat_all))
  if (length(present_targets) == 0) {
    example_ids <- head(rownames(path_mat_all), 30)
    stop(
      "None of the target pathways were found in the unstratified table rownames.\n",
      "First 30 pathway IDs seen are:\n  ",
      paste(example_ids, collapse = ", ")
    )
  }
}
message("Target pathways found (n=", length(present_targets), "): ", paste(present_targets, collapse = ", "))

## ------------------ BOX PLOTS (selected pathways; Inu only) ------------------
box_long <- as.data.frame(path_mat_all[present_targets, meta_inu$sample_id, drop = FALSE]) %>%
  rownames_to_column("pathway_id") %>%
  pivot_longer(cols = -pathway_id, names_to = "sample_id", values_to = "abundance") %>%
  left_join(meta_inu, by = "sample_id") %>%
  mutate(
    group = factor(group, levels = inu_groups),
    log10_abund = log10(abundance + 1),
    pathway_id = factor(pathway_id, levels = present_targets)
  )

box_pdf <- file.path(out_dir, "Part2_boxplots_selectedPathways_InuOnly.pdf")
safe_pdf_device(box_pdf, width = 12, height = 9)
p_box <- ggplot(box_long, aes(x = group, y = log10_abund, fill = group)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.65) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 0.8) +
  facet_wrap(~ pathway_id, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Selected pathways (Inu groups only)",
       x = "Group", y = "log10(CPM + 1)")
print(p_box)
dev.off()
message("Saved: ", box_pdf)

## ------------------ FOREST-STYLE (mean + mean-diff with bootstrap CI) ------------------
means_df <- box_long %>%
  group_by(pathway_id, group) %>%
  summarise(mean_cpm = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(panel = "Mean (CPM)", y = as.character(group), x = mean_cpm)

cmp_list <- list(
  "Inu_LP - Inu_Ctrl"   = c("Inu_Ctrl", "Inu_LP"),
  "Inu_LPI - Inu_LP"    = c("Inu_LP", "Inu_LPI"),
  "Inu_LPI - Inu_Ctrl"  = c("Inu_Ctrl", "Inu_LPI")
)

diff_rows <- list()
for (pw in present_targets) {
  df_pw <- box_long %>% filter(pathway_id == pw)
  for (cmp_name in names(cmp_list)) {
    g1 <- cmp_list[[cmp_name]][1]
    g2 <- cmp_list[[cmp_name]][2]
    x1 <- df_pw$abundance[df_pw$group == g1]
    x2 <- df_pw$abundance[df_pw$group == g2]
    ci <- boot_diff_ci(x1, x2, B = 2000, seed = 123)
    diff_rows[[length(diff_rows) + 1]] <- data.frame(
      pathway_id = pw,
      comparison = cmp_name,
      diff = as.numeric(ci["diff"]),
      lo   = as.numeric(ci["lo"]),
      hi   = as.numeric(ci["hi"]),
      panel = "Mean difference (CPM)"
    )
  }
}
diff_df <- bind_rows(diff_rows) %>%
  mutate(y = factor(comparison, levels = rev(names(cmp_list))), x = diff)

forest_pdf <- file.path(out_dir, "Part2_forestStyle_selectedPathways_InuOnly.pdf")
safe_pdf_device(forest_pdf, width = 14, height = 10)

p_forest <- ggplot() +
  geom_col(data = means_df, aes(x = x, y = y, fill = group),
           width = 0.7, alpha = 0.8) +
  geom_vline(data = diff_df %>% distinct(pathway_id, panel),
             aes(xintercept = 0),
             linewidth = 0.4, linetype = "dashed", alpha = 0.6) +
  geom_errorbarh(data = diff_df, aes(xmin = lo, xmax = hi, y = y),
                 height = 0.2, linewidth = 0.5, alpha = 0.9) +
  geom_point(data = diff_df, aes(x = x, y = y),
             size = 1.8) +
  facet_grid(pathway_id ~ panel, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 9),
        legend.position = "bottom") +
  labs(title = "Selected pathways: mean by group + mean-differences (bootstrap 95% CI)",
       x = NULL, y = NULL, fill = "Group")

print(p_forest)
dev.off()
message("Saved: ", forest_pdf)

## ------------------ STRATIFIED -> Sankey data ------------------
strat_df <- readr::read_tsv(stratified_file, show_col_types = FALSE)
strat_col <- colnames(strat_df)[1]

strat_clean <- strat_df %>%
  filter(!str_detect(.data[[strat_col]], "UNMAPPED|UNINTEGRATED|UNGROUPED"))

strat_long <- strat_clean %>%
  separate(col = !!sym(strat_col),
           into = c("pathway_full", "taxon"),
           sep = "\\|", fill = "right") %>%
  filter(!is.na(taxon)) %>%
  mutate(pathway_id = normalize_pathway_id(pathway_full)) %>%
  pivot_longer(cols = -c(pathway_full, pathway_id, taxon),
               names_to = "sample_id", values_to = "abundance") %>%
  filter(sample_id %in% meta_inu$sample_id) %>%
  left_join(meta_inu %>% select(sample_id, group), by = "sample_id") %>%
  filter(group %in% inu_groups) %>%
  mutate(genus = extract_genus(taxon)) %>%
  filter(pathway_id %in% present_targets) %>%
  filter(is.finite(abundance), abundance >= 0) %>%
  filter(!is_helicobacter(genus))

pathway_global_totals <- strat_long %>%
  group_by(pathway_id) %>%
  summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total))

pathway_levels <- pathway_global_totals$pathway_id
pathway_colors <- setNames(distinct_cols(length(pathway_levels)), pathway_levels)

genus_global_levels <- strat_long %>%
  group_by(genus) %>%
  summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genus)

prep_sankey <- function(df_long, group_name, top_genus = 25, min_flow = 0, remove_uncls = FALSE) {
  df_g <- df_long %>% filter(group == group_name)

  if (remove_uncls) {
    df_g <- df_g %>%
      filter(!is_unclassified(genus)) %>%
      filter(!is_unclassified(taxon))
  }

  agg <- df_g %>%
    group_by(pathway_id, genus) %>%
    summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(total_abundance), total_abundance > min_flow)

  if (nrow(agg) == 0) return(NULL)

  topg <- agg %>%
    group_by(genus) %>%
    summarise(gsum = sum(total_abundance), .groups = "drop") %>%
    arrange(desc(gsum)) %>%
    slice_head(n = top_genus) %>%
    pull(genus)

  agg %>%
    filter(genus %in% topg) %>%
    mutate(
      group = group_name,
      axis1 = factor(pathway_id, levels = pathway_levels),
      axis2 = factor(genus, levels = genus_global_levels),
      log_flow = log10(total_abundance + 1)
    ) %>%
    arrange(axis1, axis2)
}

plot_sankey <- function(sank_df, title, show_legend = TRUE) {
  if (is.null(sank_df) || nrow(sank_df) == 0) return(NULL)

  ggplot(sank_df, aes(axis1 = axis1, axis2 = axis2, y = log_flow)) +
    ggalluvial::geom_alluvium(aes(fill = axis1),
                              alpha = 0.80,
                              width = 1/14) +
    ggalluvial::geom_stratum(alpha = 0.95, width = 1/10) +
    ggplot2::geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 2.4
    ) +
    scale_x_discrete(limits = c("Pathway", "Genus"), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = pathway_colors, drop = FALSE) +
    theme_bw() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.title = element_blank()
    ) +
    labs(
      title = title,
      y = "log10(total stratified CPM + 1)",
      fill = "Pathway"
    )
}

make_sankey_set <- function(remove_uncls = FALSE, suffix = "") {
  ## Per-group PDFs
  for (g in inu_groups) {
    sank_g <- prep_sankey(strat_long, g, top_genus = 25, min_flow = 0, remove_uncls = remove_uncls)
    if (is.null(sank_g)) next

    out_pdf <- file.path(out_dir, paste0("Part2_sankey_pathway_to_genus_", g, suffix, ".pdf"))
    safe_pdf_device(out_pdf, width = 13, height = 9)
    print(plot_sankey(
      sank_g,
      paste0(
        "Sankey (", g, ")",
        if (remove_uncls) " - NO Unclassified" else " - WITH Unclassified",
        " | Helicobacter removed"
      ),
      show_legend = TRUE
    ))
    dev.off()
    message("Saved: ", out_pdf)
  }

  ## Multipanel comparison PDF
  sank_all <- bind_rows(lapply(inu_groups, function(g) {
    prep_sankey(strat_long, g, top_genus = 25, min_flow = 0, remove_uncls = remove_uncls)
  }))

  if (!is.null(sank_all) && nrow(sank_all) > 0) {
    multi_pdf <- file.path(out_dir, paste0("Part2_sankey_pathway_to_genus_InuGroups_MULTIPANEL", suffix, ".pdf"))
    safe_pdf_device(multi_pdf, width = 14, height = 10)

    p_multi <- ggplot(sank_all, aes(axis1 = axis1, axis2 = axis2, y = log_flow)) +
      ggalluvial::geom_alluvium(aes(fill = axis1), alpha = 0.80, width = 1/14) +
      ggalluvial::geom_stratum(alpha = 0.95, width = 1/10) +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.2) +
      facet_wrap(~ group, ncol = 1, scales = "free_y") +
      scale_x_discrete(limits = c("Pathway", "Genus"), expand = c(0.08, 0.08)) +
      scale_fill_manual(values = pathway_colors, drop = FALSE) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.title = element_blank(),
        strip.text = element_text(size = 11)
      ) +
      labs(
        title = paste0(
          "Pathway -> Genus Sankey comparison (thickness = log10(CPM+1); Helicobacter removed)",
          if (remove_uncls) " [NO Unclassified]" else " [WITH Unclassified]"
        ),
        y = "log10(total stratified CPM + 1)",
        fill = "Pathway"
      )

    print(p_multi)
    dev.off()
    message("Saved: ", multi_pdf)
  }
}

## A) WITH Unclassified
make_sankey_set(remove_uncls = FALSE, suffix = "")

## B) NO Unclassified
make_sankey_set(remove_uncls = TRUE, suffix = "_noUnclassified")

## ------------------ PIE CHARTS: Top 10 pathways per Inu group (unstratified) ------------------
TOP_N_PIE <- 10

make_pie_topN <- function(path_mat, meta_inu, group_name, topN = 10) {
  samp <- meta_inu$sample_id[meta_inu$group == group_name]
  if (length(samp) < 1) return(NULL)

  vec <- rowMeans(path_mat[, samp, drop = FALSE], na.rm = TRUE)
  tibble(pathway_id = names(vec), mean_cpm = as.numeric(vec)) %>%
    arrange(desc(mean_cpm)) %>%
    slice_head(n = topN) %>%
    mutate(
      group = group_name,
      frac = mean_cpm / sum(mean_cpm)
    )
}

pie_df <- bind_rows(lapply(inu_groups, function(g) make_pie_topN(path_mat_all, meta_inu, g, topN = TOP_N_PIE)))
if (!is.null(pie_df) && nrow(pie_df) > 0) {
  pie_pdf <- file.path(out_dir, "Part2_pie_Top10Pathways_byGroup_InuOnly.pdf")
  safe_pdf_device(pie_pdf, width = 12, height = 7)

  p_pie <- ggplot(pie_df, aes(x = "", y = frac, fill = pathway_id)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    coord_polar(theta = "y") +
    facet_wrap(~ group, ncol = 3) +
    theme_void() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 12)
    ) +
    labs(
      title = "Top 10 pathways per Inu group (unstratified; by mean CPM)",
      fill = "Pathway"
    )

  print(p_pie)
  dev.off()
  message("Saved: ", pie_pdf)
}

## =====================================================================
##  SCFA / BUTYRATE-RELATED PATHWAYS
## =====================================================================

## 1) Detect SCFA-related pathways from the FULL unstratified first-column names
scfa_ids_detected <- detect_scfa_pathways(
  pathway_full_vec = path_df_clean$pathway_full,
  pathway_id_vec   = path_df_clean$pathway_id
)
scfa_ids_detected <- intersect(scfa_ids_detected, rownames(path_mat_all))

scfa_ids <- intersect(scfa_ids_detected, rownames(path_mat_all))
if (!is.null(target_pathways)) {
  scfa_ids <- scfa_ids
}

message("SCFA/butyrate-related pathways detected (n=", length(scfa_ids), "): ",
        ifelse(length(scfa_ids) > 0, paste(scfa_ids, collapse = ", "), "NONE"))

if (length(scfa_ids) > 0) {
  scfa_box_long <- as.data.frame(path_mat_all[scfa_ids, meta_inu$sample_id, drop = FALSE]) %>%
    rownames_to_column("pathway_id") %>%
    pivot_longer(cols = -pathway_id, names_to = "sample_id", values_to = "abundance") %>%
    left_join(meta_inu, by = "sample_id") %>%
    mutate(
      group = factor(group, levels = inu_groups),
      log10_abund = log10(abundance + 1),
      pathway_id = factor(pathway_id, levels = scfa_ids)
    )

  scfa_box_pdf <- file.path(out_dir, "Part2_SCFA_boxplots_InuOnly.pdf")
  safe_pdf_device(scfa_box_pdf, width = 12, height = 9)
  p_scfa_box <- ggplot(scfa_box_long, aes(x = group, y = log10_abund, fill = group)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.65) +
    geom_jitter(width = 0.15, alpha = 0.45, size = 0.8) +
    facet_wrap(~ pathway_id, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "SCFA/butyrate-related pathways (Inu groups only)",
         x = "Group", y = "log10(CPM + 1)")
  print(p_scfa_box)
  dev.off()
  message("Saved: ", scfa_box_pdf)


  strat_scfa <- strat_long %>%
    filter(pathway_id %in% scfa_ids)

  scfa_pathway_totals <- strat_scfa %>%
    group_by(pathway_id) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total))
  scfa_path_levels <- scfa_pathway_totals$pathway_id
  scfa_path_colors <- setNames(distinct_cols(length(scfa_path_levels)), scfa_path_levels)

  prep_sankey_scfa <- function(df_long, group_name, top_genus = 25, min_flow = 0, remove_uncls = FALSE) {
    df_g <- df_long %>% filter(group == group_name)

    if (remove_uncls) {
      df_g <- df_g %>%
        filter(!is_unclassified(genus)) %>%
        filter(!is_unclassified(taxon))
    }

    agg <- df_g %>%
      group_by(pathway_id, genus) %>%
      summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      filter(is.finite(total_abundance), total_abundance > min_flow)

    if (nrow(agg) == 0) return(NULL)

    topg <- agg %>%
      group_by(genus) %>%
      summarise(gsum = sum(total_abundance), .groups = "drop") %>%
      arrange(desc(gsum)) %>%
      slice_head(n = top_genus) %>%
      pull(genus)

    agg %>%
      filter(genus %in% topg) %>%
      mutate(
        group = group_name,
        axis1 = factor(pathway_id, levels = scfa_path_levels),
        axis2 = factor(genus, levels = genus_global_levels),
        log_flow = log10(total_abundance + 1)
      ) %>%
      arrange(axis1, axis2)
  }

  plot_sankey_scfa <- function(sank_df, title, show_legend = TRUE) {
    if (is.null(sank_df) || nrow(sank_df) == 0) return(NULL)

    ggplot(sank_df, aes(axis1 = axis1, axis2 = axis2, y = log_flow)) +
      ggalluvial::geom_alluvium(aes(fill = axis1),
                                alpha = 0.82,
                                width = 1/14) +
      ggalluvial::geom_stratum(alpha = 0.95, width = 1/10) +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.4) +
      scale_x_discrete(limits = c("Pathway", "Genus"), expand = c(0.08, 0.08)) +
      scale_fill_manual(values = scfa_path_colors, drop = FALSE) +
      theme_bw() +
      theme(
        legend.position = if (show_legend) "bottom" else "none",
        axis.title = element_blank()
      ) +
      labs(
        title = title,
        y = "log10(total stratified CPM + 1)",
        fill = "Pathway"
      )
  }

  scfa_sank_all <- bind_rows(lapply(inu_groups, function(g) {
    prep_sankey_scfa(strat_scfa, g, top_genus = 25, min_flow = 0, remove_uncls = FALSE)
  }))
  if (!is.null(scfa_sank_all) && nrow(scfa_sank_all) > 0) {
    scfa_multi_pdf <- file.path(out_dir, "Part2_SCFA_sankey_PathwayToGenus_MULTIPANEL_withUnclassified.pdf")
    safe_pdf_device(scfa_multi_pdf, width = 14, height = 10)
    p_scfa_multi <- ggplot(scfa_sank_all, aes(axis1 = axis1, axis2 = axis2, y = log_flow)) +
      ggalluvial::geom_alluvium(aes(fill = axis1), alpha = 0.82, width = 1/14) +
      ggalluvial::geom_stratum(alpha = 0.95, width = 1/10) +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.2) +
      facet_wrap(~ group, ncol = 1, scales = "free_y") +
      scale_x_discrete(limits = c("Pathway", "Genus"), expand = c(0.08, 0.08)) +
      scale_fill_manual(values = scfa_path_colors, drop = FALSE) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.title = element_blank(),
        strip.text = element_text(size = 11)
      ) +
      labs(
        title = "SCFA/butyrate pathways: Pathway -> Genus (comparison across Inu groups; thickness = log10(CPM+1))",
        y = "log10(total stratified CPM + 1)",
        fill = "Pathway"
      )
    print(p_scfa_multi)
    dev.off()
    message("Saved: ", scfa_multi_pdf)
  }

  scfa_sank_all_noU <- bind_rows(lapply(inu_groups, function(g) {
    prep_sankey_scfa(strat_scfa, g, top_genus = 25, min_flow = 0, remove_uncls = TRUE)
  }))
  if (!is.null(scfa_sank_all_noU) && nrow(scfa_sank_all_noU) > 0) {
    scfa_multi_noU_pdf <- file.path(out_dir, "Part2_SCFA_sankey_PathwayToGenus_MULTIPANEL_noUnclassified.pdf")
    safe_pdf_device(scfa_multi_noU_pdf, width = 14, height = 10)
    p_scfa_multi_noU <- ggplot(scfa_sank_all_noU, aes(axis1 = axis1, axis2 = axis2, y = log_flow)) +
      ggalluvial::geom_alluvium(aes(fill = axis1), alpha = 0.82, width = 1/14) +
      ggalluvial::geom_stratum(alpha = 0.95, width = 1/10) +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.2) +
      facet_wrap(~ group, ncol = 1, scales = "free_y") +
      scale_x_discrete(limits = c("Pathway", "Genus"), expand = c(0.08, 0.08)) +
      scale_fill_manual(values = scfa_path_colors, drop = FALSE) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.title = element_blank(),
        strip.text = element_text(size = 11)
      ) +
      labs(
        title = "SCFA/butyrate pathways: Pathway -> Genus (NO Unclassified; thickness = log10(CPM+1))",
        y = "log10(total stratified CPM + 1)",
        fill = "Pathway"
      )
    print(p_scfa_multi_noU)
    dev.off()
    message("Saved: ", scfa_multi_noU_pdf)
  }

  make_group_pathway_genus_alluvial <- function(df_long, remove_uncls = FALSE, outfile) {
    dfp <- df_long

    if (remove_uncls) {
      dfp <- dfp %>%
        filter(!is_unclassified(genus)) %>%
        filter(!is_unclassified(taxon))
    }

    agg <- dfp %>%
      group_by(group, pathway_id, genus) %>%
      summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      filter(is.finite(total_abundance), total_abundance > 0) %>%
      mutate(
        group = factor(group, levels = inu_groups),
        pathway_id = factor(pathway_id, levels = scfa_path_levels),
        genus = factor(genus, levels = genus_global_levels),
        log_flow = log10(total_abundance + 1)
      )

    if (nrow(agg) == 0) return(invisible(NULL))

    safe_pdf_device(outfile, width = 15, height = 10)

    p <- ggplot(
      agg,
      aes(axis1 = group, axis2 = pathway_id, axis3 = genus, y = log_flow)
    ) +
      ggalluvial::geom_alluvium(aes(fill = pathway_id),
                                alpha = 0.80,
                                width = 1/18) +
      ggalluvial::geom_stratum(alpha = 0.95, width = 1/12) +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.2) +
      scale_x_discrete(limits = c("Group", "Pathway", "Genus"), expand = c(0.06, 0.06)) +
      scale_fill_manual(values = scfa_path_colors, drop = FALSE) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        axis.title = element_blank()
      ) +
      labs(
        title = paste0(
          "SCFA/butyrate pathways transition: Group -> Pathway -> Genus ",
          if (remove_uncls) "[NO Unclassified]" else "[WITH Unclassified]",
          " (thickness = log10(CPM+1); Helicobacter removed)"
        ),
        y = "log10(total stratified CPM + 1)",
        fill = "Pathway"
      )

    print(p)
    dev.off()
    message("Saved: ", outfile)
  }

  make_group_pathway_genus_alluvial(
    strat_scfa,
    remove_uncls = FALSE,
    outfile = file.path(out_dir, "Part2_SCFA_alluvial_Group_Pathway_Genus_withUnclassified.pdf")
  )
  make_group_pathway_genus_alluvial(
    strat_scfa,
    remove_uncls = TRUE,
    outfile = file.path(out_dir, "Part2_SCFA_alluvial_Group_Pathway_Genus_noUnclassified.pdf")
  )

} else {
  message("SCFA module skipped: no SCFA/butyrate-related pathways detected by keyword scan.")
}

## ------------------ save workspace ------------------
save.image(file = file.path(out_dir, "4b_humann_viz_part2_UPDATED_v3_scfa-ws.RData"))
message("DONE. Outputs written to: ", out_dir)
