make_fig4 <- function(out_root) {
  folder <- file.path(out_root, "fig", "Fig4_repertoire_overview")
  source_dir <- file.path(folder, "source_data")
  metrics <- read_sd(file.path(DATA,"analysis/repertoire_chain_metrics.tsv"))
  samples <- read_sd(file.path(DATA,"metadata/samples.tsv"))
  libraries <- read_sd(file.path(DATA,"metadata/libraries.tsv"))
  specimen <- c("SampleID", "tissue")
  if (!"visit_order" %in% names(metrics)) {
    metrics <- merge(metrics, samples[, c(specimen, "visit_order"), drop = FALSE],
                     by = specimen, all.x = TRUE, sort = FALSE)
  }
  metrics <- ensure_group_metadata(metrics, libraries)
  metrics$group_public <- factor_groups(metrics$group_public)
  for (x in c("strict_clonotypes", "umi_weighted_mean_cdr3_aa_length", "eligible_umi")) {
    if (x %in% names(metrics)) metrics[[x]] <- suppressWarnings(as.numeric(metrics[[x]]))
  }
  endpoint <- metrics[metrics$tissue != "blood" & metrics$processing_status == "ok", , drop = FALSE]

  native <- read_sd(file.path(out_root,"table/source_data/Table4_tissue_totals/assembled_clonotypes_by_specimen_chain.tsv"))
  native <- ensure_group_metadata(native, libraries)
  native$group_public <- factor_groups(native$group_public)
  native$chain <- factor(native$chain, levels=c("IGH","IGK","IGL","TRA","TRB","TRD","TRG"))
  native$assembled_clonotypes <- as.numeric(native$assembled_clonotypes)
  stopifnot(!anyDuplicated(native[,c(specimen,"chain")]),
            all(is.finite(native$assembled_clonotypes)), all(native$assembled_clonotypes >= 0))
  native <- merge(native, samples[,c(specimen,"visit_order","nominal_week")],
                  by=specimen, all.x=TRUE, sort=FALSE)
  endpoint_counts <- native[native$tissue != "blood", , drop=FALSE]
  write_source(endpoint_counts, file.path(source_dir, "Fig4a_assembled_clonotype_counts.tsv"))
  if (nrow(endpoint_counts)) {
    p1 <- ggplot(endpoint_counts, aes(chain, assembled_clonotypes, colour = group_public, fill = group_public, shape = group_public)) +
      geom_boxplot(aes(group = interaction(chain, group_public)), width = 0.58, alpha = 0.13, outlier.shape = NA,
                   linewidth = 0.25, position = position_dodge(width = 0.62)) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.62, seed = 20260926), size = 0.72, stroke = 0.22, alpha = 0.72) +
      facet_wrap(~tissue, scales = "free_y", nrow = 2, ncol = 5,
                 labeller = labeller(tissue = function(x) gsub("_", " ", x, fixed = TRUE))) +
      group_colour_scale(guide = "none") + group_fill_scale(guide = "none") + group_shape_scale(guide = "none") +
      scale_y_continuous(trans = "log1p", breaks = count_axis_breaks, minor_breaks = NULL,
                         labels = count_axis_labels) +
      labs(x = "Chain", y = "Clonotypes") +
      theme_sd(base_size = 7) +
      theme(axis.title = element_text(family = FIGURE_FONT, size = 7),
            axis.text = element_text(family = FIGURE_FONT, size = 6.2),
            axis.text.x = element_text(family = FIGURE_FONT, size = 6.2, angle = 90, hjust = 1, vjust = 0.5),
            strip.text = element_text(family = FIGURE_FONT, size = 7, face = "bold", margin = margin(1.5, 1, 1.5, 1)),
            panel.spacing = grid::unit(1.5, "mm"),
            plot.margin = margin(2, 2, 2, 2), legend.position = "none")
  } else p1 <- empty_panel("Clonotype counts", "No parsed endpoint sequence metrics")

  cdr3 <- endpoint[is.finite(endpoint$umi_weighted_mean_cdr3_aa_length), , drop = FALSE]
  write_source(cdr3, file.path(source_dir, "Fig4b_cdr3_length_summary.tsv"))
  if (nrow(cdr3)) {
    cdr3$chain <- factor(cdr3$chain, levels = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))
    p2 <- ggplot(cdr3, aes(group_public, umi_weighted_mean_cdr3_aa_length,
                            colour = group_public, fill = group_public, shape = group_public)) +
      geom_boxplot(width = 0.46, alpha = 0.28, outlier.shape = NA, linewidth = 0.2) +
      geom_point(position = position_jitter(width = 0.10, height = 0, seed = 20260926),
                 size = 0.5, stroke = 0.15, alpha = 0.52) +
      facet_wrap(~chain, nrow = 1) +
      group_colour_scale(name = "Group") + group_fill_scale(guide = "none") + group_shape_scale(name = "Group") +
      guides(colour = guide_legend(ncol = 1, override.aes = list(fill = unname(GROUP_COLORS)))) +
      scale_y_continuous(limits = c(10.5, 21.5), breaks = c(12, 16, 20)) +
      labs(x = NULL, y = "Mean CDR3\nlength (aa)") +
      theme_sd(base_size = 8) +
      theme(axis.text = element_text(family = FIGURE_FONT, size = 7),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title = element_text(family = FIGURE_FONT, size = 8),
            strip.text = element_text(family = FIGURE_FONT, size = 8, face = "bold",
                                      margin = margin(1, 1, 1, 1)),
            panel.spacing = grid::unit(1.2, "mm"),
            legend.position = "top", legend.direction = "horizontal",
            legend.title = element_blank(),
            legend.text = element_text(family = FIGURE_FONT, size = 8),
            legend.key.width = grid::unit(3, "mm"),
            legend.key.height = grid::unit(3, "mm"),
            plot.margin = margin(1.5, 1.5, 1.5, 1.5))
  } else p2 <- empty_panel("CDR3 length", "No parsed CDR3 length metrics")

  blood <- native[native$tissue == "blood", , drop = FALSE]
  blood$study_week <- suppressWarnings(as.numeric(blood$nominal_week))
  stopifnot(all(is.finite(blood$study_week)),
            all(blood$study_week[blood$group_public == "untreated_control"] %in% c(0,8)),
            all(blood$study_week[blood$group_public == "IMQ_treated"] %in% c(0,4,8)))
  write_source(blood, file.path(source_dir, "Fig4c_blood_native_clonotype_counts.tsv"))
  if (nrow(blood) && "visit_order" %in% names(blood)) {
    blood$visit_order <- suppressWarnings(as.numeric(blood$visit_order))
    blood$chain <- factor(blood$chain, levels = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))
    p3 <- ggplot(blood, aes(study_week, assembled_clonotypes, group = mouse_id, colour = group_public, fill = group_public, shape = group_public)) +
      geom_line(linewidth = yw_lw(0.75), alpha = 0.54) +
      geom_point(size = 0.65, stroke = 0.18, alpha = 0.82) +
      facet_wrap(~chain, nrow = 1, scales = "free_y") +
      group_colour_scale(name = "Group") + group_fill_scale(guide = "none") + group_shape_scale(name = "Group") +
      guides(colour = guide_legend(ncol = 1, override.aes = list(fill = unname(GROUP_COLORS)))) +
      scale_x_continuous(breaks = c(0, 4, 8), limits = c(-0.8, 8.8)) +
      scale_y_continuous(trans = "log1p", breaks = count_axis_breaks,
                         labels = count_axis_labels, minor_breaks = NULL) +
      labs(x = "Time (weeks)", y = "Clonotypes") +
      theme_sd(base_size = 8) +
      theme(axis.text = element_text(family = FIGURE_FONT, size = 7),
            axis.title = element_text(family = FIGURE_FONT, size = 8),
            strip.text = element_text(family = FIGURE_FONT, size = 8, face = "bold",
                                      margin = margin(1, 1, 1, 1)),
            panel.spacing = grid::unit(1.2, "mm"),
            legend.position = "top", legend.direction = "horizontal",
            legend.title = element_blank(),
            legend.text = element_text(family = FIGURE_FONT, size = 8),
            legend.key.width = grid::unit(3, "mm"),
            legend.key.height = grid::unit(3, "mm"),
            plot.margin = margin(1.5, 1.5, 1.5, 1.5))
  } else p3 <- empty_panel("Blood clonotype counts", "No blood clonotype counts supplied")

  rank_path <- file.path(DATA,"analysis/rank_abundance_top100.tsv")
  if (file.exists(rank_path) || file.exists(paste0(rank_path, ".gz"))) {
    ranks <- read_sd(rank_path)
    ranks$rank <- suppressWarnings(as.numeric(ranks$rank))
    ranks$relative_umi_abundance <- suppressWarnings(as.numeric(ranks$relative_umi_abundance))
    ranks <- ensure_group_metadata(ranks, libraries)
    ranks$group_public <- factor_groups(ranks$group_public)
    ranks <- ranks[ranks$tissue != "blood" & is.finite(ranks$rank) & is.finite(ranks$relative_umi_abundance), , drop = FALSE]
    rank_summary <- median_iqr(ranks, c("group_public", "chain", "rank"), "relative_umi_abundance")
    write_source(ranks, file.path(source_dir, "Fig4d_rank_abundance_individual_records.tsv"))
    write_source(rank_summary, file.path(source_dir, "Fig4d_rank_abundance_summary.tsv"))
    if (nrow(rank_summary)) {
      rank_summary$chain <- factor(rank_summary$chain,
                                   levels = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"))
      p4 <- ggplot(rank_summary, aes(rank, median, colour = group_public, fill = group_public, linetype = group_public)) +
        geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.17, colour = NA, show.legend = FALSE) +
        geom_line(linewidth = yw_lw(0.75)) +
        facet_wrap(~chain, nrow = 1, scales = "free_y") +
        group_colour_scale(name = "Group") + group_fill_scale(guide = "none") +
        scale_linetype_manual(values = c(untreated_control = "solid", IMQ_treated = "dashed"),
                              breaks = GROUP_ORDER, labels = GROUP_LABELS, name = "Group") +
        scale_x_continuous(breaks = c(0, 50, 100), limits = c(0, 101)) +
        scale_y_log10(labels = function(x) vapply(x, function(v) {
          if (!is.finite(v)) return(NA_character_)
          if (v >= 1e-3) format(v, scientific = FALSE, trim = TRUE)
          else sub("e-0", "e-", format(v, scientific = TRUE, digits = 1, trim = TRUE), fixed = TRUE)
        }, character(1)), breaks = function(limits) {
          10^seq.int(floor(log10(min(limits))), ceiling(log10(max(limits))))
        }) +
        labs(x = "Clonotype rank", y = "Relative UMI\nabundance") +
        theme_sd(base_size = 8) +
        theme(axis.text = element_text(family = FIGURE_FONT, size = 7),
              axis.title = element_text(family = FIGURE_FONT, size = 8),
              strip.text = element_text(family = FIGURE_FONT, size = 8, face = "bold",
                                        margin = margin(1, 1, 1, 1)),
              panel.spacing = grid::unit(1, "mm"),
              legend.position = "top", legend.direction = "horizontal", legend.title = element_blank(),
              legend.text = element_text(family = FIGURE_FONT, size = 8),
              legend.key.width = grid::unit(3, "mm"),
              legend.key.height = grid::unit(3, "mm"),
              plot.margin = margin(1.5, 1.5, 1.5, 1.5))
    } else p4 <- empty_panel("Rank-abundance profiles", "No finite rank-abundance records")
  } else {
    write_source(data.frame(), file.path(source_dir, "Fig4d_rank_abundance_summary.tsv"))
    p4 <- empty_panel("Rank-abundance profiles", "No rank_abundance_top100.tsv supplied")
  }

  legends <- c(
    "Fig. 4a. Native MiXCR clonotype counts for each endpoint specimen-chain, faceted by endpoint tissue. Counts use distinct cloneId records assigned to a single chain across the deposited exports, including productive and nonproductive records, and sum to the corresponding Table 4 tissue-group-chain totals. A natural-log ln(1+x) axis preserves zero counts while ticks retain original count units. Boxes show medians and interquartile ranges, with whiskers extending to 1.5 times that range; all observations are overlaid. No group-inference test is performed.",
    "Fig. 4b. UMI-weighted mean CDR3 amino-acid length for each endpoint sample-chain. Endpoint tissues are pooled within each chain. Boxes show medians and interquartile ranges, with whiskers extending to 1.5 times that range; all observations are overlaid. Source data retain tissue and mouse identifiers.",
    "Fig. 4c. Native MiXCR clonotype counts over time, using the same count definition as panel a and Table 4. Lines link the same mouse within each receptor chain. The x axis shows nominal study weeks (0 and 8 for controls and 0, 4 and 8 for IMQ). A natural-log ln(1+x) axis uses original count units for ticks.",
    "Fig. 4d. Median and interquartile range of relative UMI abundance at each filtered, merged clonotype rank in endpoint samples, pooling tissues within each chain. The denominator is the total UMI count of all eligible clonotypes in that specimen and chain. The abundance axis is log10. Only specimens that reach a given rank contribute at that rank. Source data retain mouse and tissue identifiers."
  )
  export_panel(p1, folder, "Fig4a_assembled_clonotype_counts", legends[1])
  export_panel(p2, folder, "Fig4b_cdr3_length_summary", legends[2])
  export_panel(p3, folder, "Fig4c_blood_native_clonotype_counts", legends[3])
  export_panel(p4, folder, "Fig4d_rank_abundance", legends[4])
  writeLines(legends, file.path(folder, "Fig4_panel_legends.txt"), useBytes = TRUE)
}
