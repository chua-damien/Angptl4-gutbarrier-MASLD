#!/usr/bin/env Rscript

## HUMAnN pathway visualisation & statistics (Path B: NO ggforce)
## FULL SCRIPT (PCA ellipse patch + ggrepel + red/white/blue heatmaps)
## + Added: PCA PC1 vs PC2 Inu-only (with envelopes)

suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "data.table", "pheatmap", "ggrepel", "circlize")
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(ggrepel)
  library(circlize)
})

## ------------------ args ------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 4a_humann_viz.R path.tsv stratified.tsv metadata.tsv [out_dir]")
}

path_file       <- args[1]
stratified_file <- args[2]
meta_file       <- args[3]
out_dir         <- ifelse(length(args) >= 4, args[4], "humann_plots")

if (!file.exists(path_file))       stop("Unstratified file not found: ", path_file)
if (!file.exists(stratified_file)) stop("Stratified file not found: ", stratified_file)
if (!file.exists(meta_file))       stop("Metadata file not found: ", meta_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ------------------ helpers ------------------

safe_pdf_device <- function(filename, width = 7, height = 6) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename, width = width, height = height)
  } else {
    grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  }
}

safe_ggsave_pdf <- function(plot, filename, width = 10, height = 8) {
  # Use cairo_pdf if available; otherwise regular pdf
  if (capabilities("cairo")) {
    ggplot2::ggsave(filename, plot = plot,
                    device = grDevices::cairo_pdf,
                    width = width, height = height)
  } else {
    ggplot2::ggsave(filename, plot = plot,
                    device = "pdf",
                    width = width, height = height)
  }
}

distinct_cols <- function(n) {
  if ("hcl.colors" %in% ls("package:grDevices")) {
    grDevices::hcl.colors(n, palette = "Dark 3")
  } else {
    grDevices::rainbow(n)
  }
}

# ---------- heatmap color palette ----------
# row-scaled heatmaps typically live around [-2, +2]
hm_breaks <- seq(-2, 2, length.out = 101)
hm_colors <- colorRampPalette(c("blue4", "white", "red4"))(length(hm_breaks) - 1)

extract_genus <- function(taxon) {
  ifelse(grepl("g__", taxon),
         sub(".*g__([^\\.\\;\\|]+).*", "\\1", taxon),
         taxon)
}

drop_unclassified <- function(df, genus_col = "genus", taxon_col = "taxon") {
  df %>%
    filter(!str_detect(.data[[genus_col]], regex("unclassified", ignore_case = TRUE))) %>%
    filter(!str_detect(.data[[taxon_col]], regex("unclassified", ignore_case = TRUE)))
}

## ---------- PCA envelopes (PATCHED: ordered + closed) ----------
ellipse_coords <- function(x, y, level = 0.95, npoints = 200) {
  x <- as.numeric(x); y <- as.numeric(y)
  df <- data.frame(x = x, y = y)
  df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (nrow(df) < 2) return(NULL)

  mu <- colMeans(df)

  theta <- seq(0, 2*pi, length.out = npoints)
  ord <- seq_along(theta)

  if (nrow(df) == 2) {
    d <- sqrt((df$x[1]-df$x[2])^2 + (df$y[1]-df$y[2])^2)
    r <- (d/2) * 1.2
    return(data.frame(
      x = mu[1] + r*cos(theta),
      y = mu[2] + r*sin(theta),
      ord = ord
    ))
  }

  S <- stats::cov(df)
  if (any(!is.finite(S)) || det(S) <= 0) {
    dmax <- max(sqrt((df$x-mu[1])^2 + (df$y-mu[2])^2))
    r <- ifelse(is.finite(dmax) && dmax > 0, dmax*1.2, 0.1)
    return(data.frame(
      x = mu[1] + r*cos(theta),
      y = mu[2] + r*sin(theta),
      ord = ord
    ))
  }

  eig <- eigen(S)
  vals <- eig$values
  vecs <- eig$vectors
  r2 <- stats::qchisq(level, df = 2)

  circle <- cbind(cos(theta), sin(theta))
  A <- vecs %*% diag(sqrt(vals * r2), 2, 2)
  pts <- t(circle %*% t(A))
  pts <- sweep(pts, 2, mu, "+")
  data.frame(x = pts[,1], y = pts[,2], ord = ord)
}

add_group_envelopes <- function(p, df, xcol, ycol, groupcol, level = 0.95) {
  df <- df %>% dplyr::filter(is.finite(.data[[xcol]]), is.finite(.data[[ycol]]))
  df[[groupcol]] <- as.factor(df[[groupcol]])

  pieces <- split(df, df[[groupcol]])
  env_list <- lapply(names(pieces), function(g) {
    d <- pieces[[g]]
    coords <- ellipse_coords(d[[xcol]], d[[ycol]], level = level)
    if (is.null(coords)) return(NULL)
    coords[[groupcol]] <- g
    coords
  })
  env <- dplyr::bind_rows(env_list)
  if (nrow(env) == 0) return(p)

  p + geom_path(
    data = env,
    aes(x = x, y = y, color = .data[[groupcol]],
        group = .data[[groupcol]], order = ord),
    linewidth = 0.7,
    linetype = "dashed",
    alpha = 0.9,
    inherit.aes = FALSE
  )
}

## ---------- Circos/chord ----------
plot_chord_with_options <- function(mat_pg, group_name, outfile_base,
                                    show_labels = FALSE,
                                    show_legends = TRUE) {
  mat_pg <- mat_pg[rowSums(mat_pg) > 0, colSums(mat_pg) > 0, drop = FALSE]
  if (nrow(mat_pg) == 0 || ncol(mat_pg) == 0) return(invisible(NULL))

  pathways <- rownames(mat_pg)
  genus <- colnames(mat_pg)

  path_cols <- distinct_cols(length(pathways)); names(path_cols) <- pathways
  gen_cols  <- distinct_cols(length(genus));    names(gen_cols)  <- genus
  grid_cols <- c(path_cols, gen_cols)

  suffix <- if (show_labels) "_withLabels" else "_noLabels"
  suffix <- paste0(suffix, if (show_legends) "_withLegend" else "_noLegend")
  outfile <- paste0(outfile_base, suffix, ".pdf")

  safe_pdf_device(outfile, width = 10, height = 9)
  on.exit({
    try(dev.off(), silent = TRUE)
    circos.clear()
  }, add = TRUE)

  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))

  if (show_labels) {
    chordDiagram(
      x = mat_pg,
      transparency = 0.35,
      annotationTrack = c("grid"),
      preAllocateTracks = list(track.height = 0.10),
      grid.col = grid_cols
    )

    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      is_path <- sector.name %in% pathways
      cex <- if (is_path) 0.45 else 0.55
      y_pos <- if (is_path) ylim[1] + 0.15 else ylim[2] - 0.15
      circos.text(
        x = mean(xlim), y = y_pos, labels = sector.name,
        facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = cex
      )
    }, bg.border = NA)

  } else {
    chordDiagram(
      x = mat_pg,
      transparency = 0.35,
      annotationTrack = c("grid"),
      preAllocateTracks = 0,
      grid.col = grid_cols
    )
  }

  title(paste0("Pathway - genus links (", group_name, ")"))

  if (show_legends) {
    par(xpd = TRUE)
    legend("left", inset = c(-0.02, 0),
           legend = pathways, fill = unname(path_cols[pathways]),
           cex = 0.6, bty = "n", title = "Pathways")
    legend("right", inset = c(-0.02, 0),
           legend = genus, fill = unname(gen_cols[genus]),
           cex = 0.7, bty = "n", title = "Genus")
    par(xpd = FALSE)
  }

  message("  Saved: ", outfile)
}

## ------------------ metadata ------------------

meta <- readr::read_tsv(meta_file, show_col_types = FALSE)
if (!all(c("sample_id", "group") %in% colnames(meta))) {
  stop("Metadata must have columns: sample_id, group")
}

group_order <- c("2w_Ctrl","2w_LP","Inu_Ctrl","Inu_LP","Inu_LPI","Rev_Ctrl","Rev_LP","Rev_R")
meta$group <- factor(meta$group, levels = group_order)

## ------------------ unstratified pathways ------------------

message("Reading unstratified HUMAnN pathway table: ", path_file)
path_df <- readr::read_tsv(path_file, show_col_types = FALSE)
path_col <- colnames(path_df)[1]

path_df_clean <- path_df %>%
  filter(
    !str_detect(.data[[path_col]], "UNMAPPED|UNINTEGRATED|UNGROUPED"),
    !str_detect(.data[[path_col]], "\\|")
  )

sample_cols <- colnames(path_df_clean)[-1]
common_samples <- intersect(sample_cols, meta$sample_id)
if (length(common_samples) < 2) {
  stop("Fewer than 2 overlapping samples between pathway table and metadata.")
}

path_df_clean <- path_df_clean %>% select(all_of(path_col), all_of(common_samples))

meta <- meta %>% filter(sample_id %in% common_samples)
meta <- meta[match(common_samples, meta$sample_id), ]
meta$group <- factor(meta$group, levels = group_order)

path_mat <- path_df_clean %>% column_to_rownames(path_col) %>% as.matrix()
path_mat_log <- log10(path_mat + 1)

## ------------------ PCA of samples ------------------

message("Running PCA on samples...")

n_top_pca <- min(500, nrow(path_mat_log))
vars <- apply(path_mat_log, 1, var)
top_ids <- names(sort(vars, decreasing = TRUE))[1:n_top_pca]
pca_input <- t(path_mat_log[top_ids, , drop = FALSE])

pca_res <- prcomp(pca_input, scale. = TRUE)
pca_scores <- as.data.frame(pca_res$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(group = as.factor(group))

pca_label_layer <- ggrepel::geom_text_repel(
  aes(label = sample_id),
  size = 2.3,
  max.overlaps = 25,
  box.padding = 0.5,
  point.padding = 0.3,
  force = 3,
  min.segment.length = 0,
  segment.alpha = 0.35,
  show.legend = FALSE
)

# PC1 vs PC2
pca12_file <- file.path(out_dir, "PCA_PC1_PC2_samples.pdf")
safe_pdf_device(pca12_file, width = 8.2, height = 6.6)

p <- ggplot(pca_scores, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.3) +
  pca_label_layer +
  theme_minimal() +
  labs(
    title = "PCA of pathways (PC1 vs PC2)",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
  )

p <- add_group_envelopes(p, pca_scores, "PC1", "PC2", "group", level = 0.95)
print(p)
dev.off()
message("Saved: ", pca12_file)

# PC1 vs PC3
if (ncol(pca_res$x) >= 3) {
  pca13_file <- file.path(out_dir, "PCA_PC1_PC3_samples.pdf")
  safe_pdf_device(pca13_file, width = 8.2, height = 6.6)

  p <- ggplot(pca_scores, aes(PC1, PC3, color = group)) +
    geom_point(size = 2.3) +
    pca_label_layer +
    theme_minimal() +
    labs(
      title = "PCA of pathways (PC1 vs PC3)",
      x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC3 (", round(summary(pca_res)$importance[2, 3] * 100, 1), "%)")
    )

  p <- add_group_envelopes(p, pca_scores, "PC1", "PC3", "group", level = 0.95)
  print(p)
  dev.off()
  message("Saved: ", pca13_file)
}

## ------------------ PCA (Inu-only; PC1 vs PC2) ------------------

message("Running PCA on Inu-only samples (PC1 vs PC2)...")

inu_only_idx <- grepl("^Inu", as.character(pca_scores$group))
pca_scores_inu <- pca_scores[inu_only_idx, , drop = FALSE]

if (nrow(pca_scores_inu) < 3 || length(unique(pca_scores_inu$group)) < 2) {
  message("Skipping Inu-only PCA: not enough samples or groups.")
} else {

  pca_inu_file <- file.path(out_dir, "PCA_PC1_PC2_InuOnly_samples.pdf")
  safe_pdf_device(pca_inu_file, width = 8.2, height = 6.6)

  p_inu <- ggplot(pca_scores_inu, aes(PC1, PC2, color = group)) +
    geom_point(size = 2.6) +
    ggrepel::geom_text_repel(
      aes(label = sample_id),
      size = 2.4,
      max.overlaps = 30,
      box.padding = 0.6,
      point.padding = 0.35,
      force = 4,
      min.segment.length = 0,
      segment.alpha = 0.35,
      show.legend = FALSE
    ) +
    theme_minimal() +
    labs(
      title = "PCA of pathways (Inu-only samples, PC1 vs PC2)",
      x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
    )

  p_inu <- add_group_envelopes(p_inu, pca_scores_inu, "PC1", "PC2", "group", level = 0.95)
  print(p_inu)
  dev.off()

  message("Saved: ", pca_inu_file)
}

message("Running pairwise Kruskal-Wallis tests...")

comparisons <- list(
  Inu_LP_vs_Inu_Ctrl  = c("Inu_LP",  "Inu_Ctrl"),
  Inu_LPI_vs_Inu_LP   = c("Inu_LPI", "Inu_LP"),
  Inu_LPI_vs_Inu_Ctrl = c("Inu_LPI", "Inu_Ctrl")
)

kw_results_list <- list()
selected_paths <- character(0)

for (cmp_name in names(comparisons)) {
  groups_pair <- comparisons[[cmp_name]]
  message("  Comparison: ", cmp_name, " (", groups_pair[1], " vs ", groups_pair[2], ")")

  idx <- meta$group %in% groups_pair
  grp_pair <- droplevels(meta$group[idx])

  if (length(unique(grp_pair)) < 2 || sum(idx) < 3) {
    message("    Skipping: not enough samples or groups.")
    next
  }

  mat_pair <- path_mat_log[, idx, drop = FALSE]
  tests <- apply(mat_pair, 1, function(x) suppressWarnings(kruskal.test(x ~ grp_pair)))
  p_vals <- sapply(tests, function(z) unname(z$p.value))
  H_vals <- sapply(tests, function(z) unname(z$statistic))
  fdr_vals <- p.adjust(p_vals, method = "BH")

  k <- length(unique(grp_pair))
  N <- length(grp_pair)
  eta2_vals <- (H_vals - (k - 1)) / (N - 1)
  eta2_vals[eta2_vals < 0] <- 0

  mean_g1 <- rowMeans(path_mat[, meta$group == groups_pair[1], drop = FALSE])
  mean_g2 <- rowMeans(path_mat[, meta$group == groups_pair[2], drop = FALSE])

  res <- data.frame(
    pathway     = rownames(path_mat_log),
    group1      = groups_pair[1],
    group2      = groups_pair[2],
    p_value     = p_vals,
    fdr         = fdr_vals,
    eta2        = eta2_vals,
    mean_group1 = mean_g1[rownames(path_mat_log)],
    mean_group2 = mean_g2[rownames(path_mat_log)],
    diff_mean   = mean_g2[rownames(path_mat_log)] - mean_g1[rownames(path_mat_log)],
    stringsAsFactors = FALSE
  )

  res <- res[order(res$fdr), ]
  csv_file <- file.path(out_dir, paste0("KW_", cmp_name, ".csv"))
  write.csv(res, csv_file, row.names = FALSE)
  message("    Saved KW table: ", csv_file)

  top50 <- head(res$pathway, 50)
  selected_paths <- union(selected_paths, top50)
  kw_results_list[[cmp_name]] <- res
}

if (length(selected_paths) == 0) stop("No pathways selected from pairwise KW tests.")
message("Heatmap pathway union size (top50 per comparison union): ", length(selected_paths))

## ------------------ Heatmaps ------------------

message("Building heatmap for UNION of top 50 pathways from EACH comparison (expect ~150)...")

selected_paths <- intersect(selected_paths, rownames(path_mat_log))
heat_mat <- path_mat_log[selected_paths, , drop = FALSE]

ord <- order(meta$group, meta$sample_id)
heat_mat <- heat_mat[, ord, drop = FALSE]
ann_col <- data.frame(group = meta$group[ord]); rownames(ann_col) <- meta$sample_id[ord]

heat_file <- file.path(out_dir, "heatmap_union_top50each_byComparison_rowScaled.pdf")
safe_pdf_device(heat_file, width = 10, height = 10)
pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_col = ann_col,
  color = hm_colors,
  breaks = hm_breaks,
  main = paste0("UNION of top50 per comparison (n = ", length(selected_paths), ") [row-scaled]"),
  fontsize_row = 4
)
dev.off()
message("Saved: ", heat_file)

message("Building Inu-only heatmap (columns where group starts with 'Inu')...")

inu_idx <- grepl("^Inu", as.character(meta$group))
if (sum(inu_idx) < 2) {
  message("Skipping Inu-only heatmap: fewer than 2 Inu samples found.")
} else {
  meta_inu <- meta[inu_idx, , drop = FALSE]
  heat_inu <- path_mat_log[selected_paths, meta_inu$sample_id, drop = FALSE]

  ord_inu <- order(meta_inu$group, meta_inu$sample_id)
  heat_inu <- heat_inu[, ord_inu, drop = FALSE]
  ann_inu <- data.frame(group = meta_inu$group[ord_inu]); rownames(ann_inu) <- meta_inu$sample_id[ord_inu]

  heat_inu_file <- file.path(out_dir, "heatmap_union_top50each_InuOnly_rowScaled.pdf")
  safe_pdf_device(heat_inu_file, width = 9, height = 10)
  pheatmap(
    heat_inu,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    annotation_col = ann_inu,
    color = hm_colors,
    breaks = hm_breaks,
    main = paste0("UNION top50 per comparison (n = ", length(selected_paths), ") | Inu-only [row-scaled]"),
    fontsize_row = 4
  )
  dev.off()
  message("Saved: ", heat_inu_file)
}

## ------------------ Forest plots ------------------

message("Building effect-size vs FDR plots per comparison (repelled labels + zoomed-out version)...")

N_LABELS_FOCUSED <- 25
N_LABELS_ZOOMED  <- 40

for (cmp_name in names(kw_results_list)) {
  res <- kw_results_list[[cmp_name]]
  res <- res[is.finite(res$fdr) & !is.na(res$fdr), ]
  res$neglog <- -log10(res$fdr + 1e-300)

  res <- res %>% arrange(fdr, desc(eta2))

  # Zoomed labels
  res$label <- ""
  nZ <- min(N_LABELS_ZOOMED, nrow(res))
  if (nZ > 0) res$label[seq_len(nZ)] <- res$pathway[seq_len(nZ)]

  # Focused labels
  res_focus <- res
  res_focus$label <- ""
  nF <- min(N_LABELS_FOCUSED, nrow(res_focus))
  if (nF > 0) res_focus$label[seq_len(nF)] <- res_focus$pathway[seq_len(nF)]

  forest_file <- file.path(out_dir, paste0("forest_", cmp_name, "_focused.pdf"))
  safe_pdf_device(forest_file, width = 9.5, height = 7.0)

  p1 <- ggplot(res_focus, aes(x = neglog, y = eta2)) +
    geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 2.4,
      max.overlaps = Inf,
      box.padding = 0.7,
      point.padding = 0.4,
      force = 6,
      force_pull = 0.5,
      min.segment.length = 0,
      segment.alpha = 0.35,
      show.legend = FALSE
    ) +
    theme_minimal() +
    labs(
      x = "-log10(FDR)",
      y = "Effect size (eta2)",
      title = paste0("Pathway effect size vs significance: ", cmp_name, " (focused)")
    )
  print(p1)
  dev.off()
  message("Saved: ", forest_file)

  forest_file2 <- file.path(out_dir, paste0("forest_", cmp_name, "_zoomedOut.pdf"))
  safe_pdf_device(forest_file2, width = 11.0, height = 8.5)

  x_max <- max(res$neglog, na.rm = TRUE)
  y_max <- max(res$eta2, na.rm = TRUE)

  p2 <- ggplot(res, aes(x = neglog, y = eta2)) +
    geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 2.2,
      max.overlaps = Inf,
      box.padding = 0.9,
      point.padding = 0.45,
      force = 6,
      force_pull = 0.5,
      min.segment.length = 0,
      segment.alpha = 0.3,
      show.legend = FALSE
    ) +
    coord_cartesian(
      xlim = c(0, x_max * 1.15),
      ylim = c(0, y_max * 1.15)
    ) +
    theme_minimal() +
    labs(
      x = "-log10(FDR)",
      y = "Effect size (eta2)",
      title = paste0("Pathway effect size vs significance: ", cmp_name, " (zoomed-out)")
    )
  print(p2)
  dev.off()
  message("Saved: ", forest_file2)
}

## ------------------ Stratified: chord plots + boxplots ------------------

message("Reading stratified HUMAnN table: ", stratified_file)

strat_df <- readr::read_tsv(stratified_file, show_col_types = FALSE)
strat_col <- colnames(strat_df)[1]

strat_clean <- strat_df %>%
  filter(!str_detect(.data[[strat_col]], "UNMAPPED|UNINTEGRATED|UNGROUPED"))

strat_long <- strat_clean %>%
  separate(col = !!sym(strat_col),
           into = c("pathway", "taxon"),
           sep = "\\|",
           fill = "right") %>%
  filter(!is.na(taxon)) %>%
  pivot_longer(
    cols = -c(pathway, taxon),
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  filter(sample_id %in% meta$sample_id)

if (nrow(strat_long) == 0) {
  message("No stratified rows after filtering; skipping chord/circos and boxplots.")
} else {

  strat_long <- strat_long %>% mutate(genus = extract_genus(taxon))
  target_groups <- c("Inu_Ctrl", "Inu_LP", "Inu_LPI")

  for (g in target_groups) {
    message("Building chord plots for group: ", g)

    samples_g <- meta$sample_id[meta$group == g]
    if (length(samples_g) == 0) { message("  No samples for group ", g, "; skipping."); next }

    df_g <- strat_long %>% filter(sample_id %in% samples_g)
    if (nrow(df_g) == 0) { message("  No stratified entries for group ", g, "; skipping."); next }

    ps_agg <- df_g %>%
      group_by(pathway, genus) %>%
      summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      filter(is.finite(total_abundance), total_abundance > 0)

    if (nrow(ps_agg) == 0) { message("  No valid abundance entries for group ", g, "; skipping."); next }

    top_pathways <- ps_agg %>%
      group_by(pathway) %>%
      summarise(psum = sum(total_abundance), .groups = "drop") %>%
      arrange(desc(psum)) %>%
      slice_head(n = 20) %>%
      pull(pathway)

    top_genus <- ps_agg %>%
      group_by(genus) %>%
      summarise(gsum = sum(total_abundance), .groups = "drop") %>%
      arrange(desc(gsum)) %>%
      slice_head(n = 15) %>%
      pull(genus)

    chord_df <- ps_agg %>% filter(pathway %in% top_pathways, genus %in% top_genus)
    if (nrow(chord_df) > 0) {
      mat_pg <- as.matrix(xtabs(total_abundance ~ pathway + genus, data = chord_df))
      outfile_base <- file.path(out_dir, paste0("chord_pathway_genus_", g))

      plot_chord_with_options(mat_pg, g, outfile_base, show_labels = FALSE, show_legends = TRUE)
      plot_chord_with_options(mat_pg, g, outfile_base, show_labels = TRUE,  show_legends = FALSE)
      plot_chord_with_options(mat_pg, g, outfile_base, show_labels = TRUE,  show_legends = TRUE)
    }

    # repeat after removing Unclassified
    df_g_noU <- df_g %>% drop_unclassified(genus_col = "genus", taxon_col = "taxon")
    ps_agg2 <- df_g_noU %>%
      group_by(pathway, genus) %>%
      summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      filter(is.finite(total_abundance), total_abundance > 0)

    if (nrow(ps_agg2) == 0) { message("  After removing Unclassified, no data for group ", g, "; skipping."); next }

    top_genus2 <- ps_agg2 %>%
      filter(pathway %in% top_pathways) %>%
      group_by(genus) %>%
      summarise(gsum = sum(total_abundance), .groups = "drop") %>%
      arrange(desc(gsum)) %>%
      slice_head(n = 15) %>%
      pull(genus)

    chord_df2 <- ps_agg2 %>% filter(pathway %in% top_pathways, genus %in% top_genus2)
    if (nrow(chord_df2) > 0) {
      mat_pg2 <- as.matrix(xtabs(total_abundance ~ pathway + genus, data = chord_df2))
      outfile_base2 <- file.path(out_dir, paste0("chord_pathway_genus_", g, "_noUnclassified"))

      plot_chord_with_options(mat_pg2, paste0(g, ", no Unclassified"), outfile_base2, show_labels = FALSE, show_legends = TRUE)
      plot_chord_with_options(mat_pg2, paste0(g, ", no Unclassified"), outfile_base2, show_labels = TRUE,  show_legends = FALSE)
      plot_chord_with_options(mat_pg2, paste0(g, ", no Unclassified"), outfile_base2, show_labels = TRUE,  show_legends = TRUE)
    }
  }

  ## ---- boxplots of top taxa by group (stratified) ----
  message("Building boxplots of top taxa by group...")

  strat_sample <- strat_long %>%
    group_by(taxon, sample_id) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    left_join(meta, by = "sample_id") %>%
    filter(is.finite(abundance), abundance >= 0)

  if (nrow(strat_sample) == 0) {
    message("No valid stratified taxon abundance values; skipping boxplots.")
  } else {
    top_species <- strat_sample %>%
      group_by(taxon) %>%
      summarise(tot = sum(abundance), .groups = "drop") %>%
      arrange(desc(tot)) %>%
      slice_head(n = 9) %>%
      pull(taxon)

    box_df <- strat_sample %>%
      filter(taxon %in% top_species) %>%
      mutate(
        taxon = factor(taxon, levels = top_species),
        log_abundance = log10(abundance + 1),
        group = factor(group, levels = group_order)
      ) %>%
      filter(is.finite(log_abundance))

    if (nrow(box_df) == 0) {
      message("Boxplot dataframe empty after filtering; skipping boxplots.")
    } else {
      p_box <- ggplot(box_df, aes(x = group, y = log_abundance, fill = group)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.6) +
        geom_jitter(width = 0.15, alpha = 0.4, size = 0.6) +
        facet_wrap(~ taxon, scales = "free_y") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(x = "Group", y = "log10(relative abundance + 1)",
             title = "Top taxa abundance across groups")

      box_file <- file.path(out_dir, "boxplots_topTaxa_byGroup.pdf")
      safe_ggsave_pdf(p_box, box_file, width = 11, height = 8)
      message("Saved: ", box_file)
    }
  }
}

save.image(file = file.path(out_dir, "4a_humann_viz-ws.RData"))

message("All plots written to: ", out_dir)
