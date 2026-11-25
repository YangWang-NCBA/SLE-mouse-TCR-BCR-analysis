make_fig2 <- function(out_root) {
out <- file.path(out_root, "fig", "Fig2_sequencing_quality")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
src <- file.path(out, "source_data")
dir.create(src,recursive=TRUE,showWarnings=FALSE)
read_table <- function(name) {
  generated <- file.path(src,name)
  read_sd(if (file.exists(generated)) generated else file.path(DATA,"fastqc",name))
}
write_table <- function(x,name) write_source(x,file.path(src,name))
groups <- GROUP_ORDER
colors <- GROUP_COLORS
labels <- GROUP_LABELS
ink <- "#263238"
lw <- function(pt) pt / (0.75 * ggplot2::.pt)

group_color <- function() scale_color_manual(values = colors, breaks = groups, labels = labels, drop = FALSE)
group_fill <- function() scale_fill_manual(values = colors, breaks = groups, labels = labels, drop = FALSE)
group_shape <- function() scale_shape_manual(values = c(untreated_control = 1, IMQ_treated = 2), breaks = groups, labels = labels)
summarize_bins <- function(data, value, keys) {
  spl <- split(data, interaction(data[keys], drop = TRUE, lex.order = TRUE))
  ans <- lapply(spl, function(x) {
    v <- x[[value]]
    stopifnot(all(is.finite(v)))
    cbind(x[1, keys, drop = FALSE], data.frame(n = length(v), median = median(v), q25 = unname(quantile(v, 0.25)), q75 = unname(quantile(v, 0.75))))
  })
  ans <- do.call(rbind, ans); rownames(ans) <- NULL
  ans[do.call(order, ans[keys]), , drop = FALSE]
}
per_base <- read_table("raw_per_base_metrics.tsv")
q <- summarize_bins(per_base, "mean_quality", c("group_public", "read_direction", "base_midpoint"))
write_table(q,"Fig2a_fastqc_per_base_quality_summary.tsv")
q$group_public <- factor(q$group_public, groups)
p1 <- ggplot(q, aes(base_midpoint, median, color = group_public, fill = group_public,
                   linetype = read_direction, group = interaction(group_public, read_direction))) +
  geom_ribbon(aes(ymin = q25, ymax = q75, linetype = NULL), alpha = 0.13, color = NA) +
  geom_line(linewidth = lw(0.75)) +
  geom_hline(yintercept = 30, color = "#7A7A7A", linetype = "dotted", linewidth = lw(0.5)) +
  group_color() + group_fill() + scale_linetype_manual(values = c(R1 = "solid", R2 = "22")) +
  scale_y_continuous(breaks = seq(20, 40, 5), expand = expansion(mult = 0)) +
  coord_cartesian(ylim = c(20, 40)) +
  guides(fill = "none", color = guide_legend(order = 1), linetype = guide_legend(order = 2, override.aes = list(fill = NA))) +
  labs(x = "Read position (bp)", y = "Mean Phred quality", tag = "a") + theme_sd(base_size = 6, compact = TRUE)
pair <- read_table("Fig2b_R1_R2_read_count_agreement.tsv")
pair$group_public <- factor(pair$group_public, groups)
p2 <- ggplot(pair, aes(R1_read_count, R2_read_count, color = group_public, shape = group_public)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#7A7A7A", linewidth = lw(0.5)) +
  geom_point(size = 1.15, stroke = 0.35, alpha = 0.76) + group_color() + group_shape() +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "R1 reads", y = "R2 reads", tag = "b") + theme_sd(base_size = 6, compact = TRUE)
n <- per_base
stopifnot(all(n$n_percent >= 0 & n$n_percent <= 100))
n$group_public <- factor(n$group_public, groups)
n <- n[order(n$filename, n$base_midpoint), ]
ns <- summarize_bins(n, "n_percent", c("group_public", "read_direction", "base_midpoint"))
write_table(ns, "Fig2c_per_base_n_content_summary.tsv")
p3 <- ggplot() +
  geom_line(data = n, aes(base_midpoint, n_percent, group = filename, color = group_public), linewidth = lw(0.2), alpha = 0.065) +
  geom_ribbon(data = ns, aes(base_midpoint, ymin = q25, ymax = q75, fill = group_public, group = group_public), alpha = 0.16) +
  geom_line(data = ns, aes(base_midpoint, median, color = group_public, group = group_public), linewidth = lw(0.9)) +
  facet_wrap(~read_direction, nrow = 1) + group_color() + group_fill() +
  scale_x_continuous(breaks = c(0, 150, 300), limits = c(0, 300), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, ceiling(max(n$n_percent))), expand = expansion(mult = c(0, 0.035))) +
  guides(fill = "none", color = guide_legend(override.aes = list(alpha = 1, linewidth = lw(0.9)))) +
  labs(x = "Read position (bp)", y = "N content (%)", tag = "c") + theme_sd(base_size = 6, compact = TRUE)

d <- read_table("Fig2d_trimmed_read_length_distribution.tsv")
# Exact 1-bp counts are normalized within each file. Unobserved lengths
# were filled with zero before the equal-library group average is computed.
stopifnot(all(c("analysis_library_id", "read_direction", "receptor", "group_public", "length_midpoint", "percent_per_bp") %in% names(d)))
stopifnot(all(abs(aggregate(percent_per_bp ~ analysis_library_id + read_direction, d, sum)$percent_per_bp - 100) < 1e-8))
d$group_public <- factor(d$group_public, groups)
ds <- aggregate(percent_per_bp ~ receptor + read_direction + group_public + length_midpoint, d, mean)
write_table(ds, "Fig2d_trimmed_read_length_summary.tsv")
p4 <- ggplot(ds, aes(length_midpoint, percent_per_bp, color = group_public, linetype = group_public)) +
  geom_step(linewidth = lw(0.75), direction = "mid") +
  facet_grid(receptor ~ read_direction) + group_color() +
  scale_linetype_manual(values = c(untreated_control = "solid", IMQ_treated = "22"), breaks = groups, labels = labels) +
  scale_x_continuous(breaks = c(0, 150, 300), limits = c(0, 300)) +
  labs(x = "Read length after trimming (bp)", y = "Reads (%)", tag = "d") + theme_sd(base_size = 6, compact = TRUE)

draw_four <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, heights = grid::unit(c(0.49, 0.51), "null"))))
  for (i in seq_len(4)) print(list(p1, p2, p3, p4)[[i]], newpage = FALSE,
    vp = grid::viewport(layout.pos.row = if (i <= 2) 1 else 2, layout.pos.col = if (i %% 2) 1 else 2))
  grid::popViewport()
}
stem <- "Fig2_assembled"
cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = 183 / 25.4, height = 150 / 25.4, family = FIGURE_FONT, bg = "white")
draw_four(); invisible(dev.off())
if ("svg" %in% FORMATS) {
svglite::svglite(file.path(out, paste0(stem, ".svg")), width = 183 / 25.4, height = 150 / 25.4, bg = "white", system_fonts = list(sans = FIGURE_FONT))
draw_four(); invisible(dev.off())
}
png(file.path(out, paste0(stem, ".png")), width = 183, height = 150, units = "mm", res = 600, type = "cairo", bg = "white")
draw_four(); invisible(dev.off())
panel_notes <- c("Fig. 2a | Raw-read per-base quality. At each position, lines show the median per-file mean Phred score and shading shows the interquartile range within each treatment group and read direction (120 control and 130 IMQ files per direction). Solid and dashed lines denote R1 and R2, respectively; the dotted line marks Q30.",
  "Fig. 2b | Raw R1/R2 read-count agreement. Each point is one of 250 paired-end libraries. R1 and R2 counts agree in every library; the diagonal denotes equal counts.",
  "Fig. 2c | Raw-read per-base N content. Thin lines represent all 500 files, with median and interquartile range by group and read direction. No files are excluded.",
  "Fig. 2d | Read lengths after paired-end adapter and Q20 quality trimming, with a minimum length of 35 bp. Counts in 1-bp bins are normalized within each file and averaged with equal library weights, separately by group, receptor type and read direction.")
for (i in seq_len(4)) {
  spec_row <- PANEL_EXPORT_SPEC[PANEL_EXPORT_SPEC$panel_id == paste0("Fig2", letters[i]), ]
  export_panel(list(p1,p2,p3,p4)[[i]], out, spec_row$stem, panel_notes[i])
}

invisible(out)
}
