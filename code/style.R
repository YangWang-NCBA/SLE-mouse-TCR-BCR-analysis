FIGURE_FONT <- Sys.getenv("SD_FONT_FAMILY", unset = "Arial")

xml_font_path <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

configure_figure_fonts <- function() {
  if (!isTRUE(capabilities("cairo"))) stop("R requires Cairo support for PDF/PNG/TIFF export", call. = FALSE)
  font_dir <- if (FIGURE_FONT == "Arial") Sys.getenv("ARIAL_FONT_DIR", unset = "") else ""
  windows <- identical(.Platform$OS.type, "windows")
  if (nzchar(font_dir)) {
    font_dir <- normalizePath(font_dir, winslash = "/", mustWork = TRUE)
    if (!dir.exists(font_dir) || !length(list.files(font_dir, pattern = "[.](ttf|otf)$", ignore.case = TRUE))) {
      stop("ARIAL_FONT_DIR must contain licensed Arial font files", call. = FALSE)
    }
    if (!windows) {
      cache_dir <- file.path(tempdir(), "fontconfig_cache")
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      config_path <- file.path(tempdir(), "fontconfig.xml")
      writeLines(c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">',
        '<fontconfig>', '<include ignore_missing="yes">/etc/fonts/fonts.conf</include>',
        paste0('<dir>', xml_font_path(font_dir), '</dir>'),
        paste0('<cachedir>', xml_font_path(normalizePath(cache_dir, winslash = "/", mustWork = TRUE)), '</cachedir>'),
        '</fontconfig>'), config_path, useBytes = TRUE)
      Sys.setenv(FONTCONFIG_FILE = normalizePath(config_path, winslash = "/", mustWork = TRUE))
    }
  }
  handoff <- Sys.getenv("SD_BOLD_FONT_FILE")
  if (nzchar(handoff)) {
    face <- systemfonts::match_fonts(FIGURE_FONT, weight = "bold")
    info <- systemfonts::font_info(path = face$path[1], index = face$index[1])
    if (!identical(tolower(as.character(info$family[1])), tolower(FIGURE_FONT)) || !isTRUE(as.logical(info$bold[1]))) {
      stop("Install the selected bold font: ", FIGURE_FONT, call. = FALSE)
    }
    writeLines(enc2utf8(face$path[1]), handoff, useBytes = TRUE)
  }
}

theme_sd <- function(base_size = 10, compact = FALSE) {
  base <- if (compact) theme_classic(base_family = FIGURE_FONT, base_size = base_size) else
    theme_minimal(base_family = FIGURE_FONT, base_size = base_size)
  result <- base + theme(
    text = element_text(family = FIGURE_FONT, colour = INK),
    axis.title = element_text(family = FIGURE_FONT, size = if (compact) base_size + 1 else base_size,
                             margin = if (compact) NULL else margin(t = 5, r = 5)),
    axis.text = element_text(family = FIGURE_FONT, size = if (compact) base_size else 9.5, colour = INK),
    strip.text = element_text(family = FIGURE_FONT, size = if (compact) base_size else 10,
                             face = if (compact) "plain" else "bold",
                             margin = if (compact) margin(3, 2, 3, 2) else margin(4, 4, 4, 4)),
    strip.background = element_rect(fill = "#F2F4F4", colour = NA),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.spacing = grid::unit(if (compact) 3.5 else 5, "mm"),
    legend.title = if (compact) element_blank() else element_text(family = FIGURE_FONT, size = 9.5),
    legend.text = element_text(family = FIGURE_FONT, size = if (compact) base_size else 9.5),
    legend.position = "bottom", legend.box = if (compact) "horizontal" else "vertical",
    legend.spacing.x = grid::unit(if (compact) 1.5 else 2.5, "mm"),
    legend.key.width = grid::unit(if (compact) 5 else 6, "mm"),
    legend.key.height = grid::unit(if (compact) 3.2 else 4, "mm"),
    plot.tag = element_blank(),
    plot.margin = if (compact) margin(5, 5, 4, 5) else margin(8, 12, 8, 8)
  )
  if (compact) {
    result <- result + theme(
      text = element_text(face = "plain"), axis.title = element_text(face = "plain"),
      axis.line = element_line(colour = "black", linewidth = yw_lw(0.5)),
      axis.ticks = element_line(colour = "black", linewidth = yw_lw(0.5)),
      axis.ticks.length = grid::unit(2.5, "pt"),
      panel.grid.major.y = element_line(colour = "#E7EBEF", linewidth = yw_lw(0.3)),
      legend.margin = margin(1, 0, 0, 0),
      plot.tag = element_text(size = 8, face = "bold", family = FIGURE_FONT, colour = "black"),
      plot.tag.position = c(0, 1), plot.tag.location = "plot"
    )
  } else {
    result <- result + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                             plot.caption = element_blank())
  }
  result
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

GROUP_ORDER <- c("untreated_control", "IMQ_treated")
GROUP_COLORS <- c(untreated_control = "#438BB5", IMQ_treated = "#E88F86")
GROUP_SHAPES <- c(untreated_control = 21, IMQ_treated = 24)
GROUP_LABELS <- c(untreated_control = "Control", IMQ_treated = "IMQ")
INK <- "#263238"
MID <- "#7A7A7A"

read_sd <- function(path) {
  if (!file.exists(path) && file.exists(paste0(path, ".gz"))) path <- paste0(path, ".gz")
  if (!file.exists(path)) stop("Required input missing: ", path, call. = FALSE)
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", ""), fileEncoding = "UTF-8")
}

factor_groups <- function(x) factor(x, levels = GROUP_ORDER)

ensure_group_metadata <- function(data, metadata = NULL, key = "analysis_library_id") {
  groups <- if ("group_public" %in% names(data)) as.character(data$group_public) else rep(NA_character_, nrow(data))
  if (!is.null(metadata)) {
    if (!key %in% names(data) || !all(c(key, "group_public") %in% names(metadata))) {
      stop("Group metadata requires ", key, " in the data and metadata, and group_public in the metadata", call. = FALSE)
    }
    lookup <- unique(metadata[, c(key, "group_public"), drop = FALSE])
    lookup[[key]] <- as.character(lookup[[key]])
    lookup$group_public <- as.character(lookup$group_public)
    if (any(is.na(lookup[[key]]) | !nzchar(lookup[[key]])) ||
        any(is.na(lookup$group_public) | !lookup$group_public %in% GROUP_ORDER)) {
      stop("Group metadata contains missing identifiers or unsupported group labels", call. = FALSE)
    }
    if (anyDuplicated(lookup[[key]])) {
      stop("Group metadata has conflicting group labels for ", key, call. = FALSE)
    }
    index <- match(as.character(data[[key]]), lookup[[key]])
    if (anyNA(index)) {
      stop("Plotting input contains identifiers absent from the group metadata: ", key, call. = FALSE)
    }
    expected <- lookup$group_public[index]
    supplied <- !is.na(groups) & nzchar(groups)
    if (any(supplied & groups != expected)) {
      stop("Plotting group_public disagrees with the explicit group metadata", call. = FALSE)
    }
    groups[!supplied] <- expected[!supplied]
  }
  if (any(is.na(groups) | !groups %in% GROUP_ORDER)) {
    stop("Plotting requires valid group_public values or explicit group metadata", call. = FALSE)
  }
  data$group_public <- groups
  data
}

group_colour_scale <- function(...) scale_colour_manual(values = GROUP_COLORS, breaks = GROUP_ORDER, labels = GROUP_LABELS, drop = FALSE, ...)
group_fill_scale <- function(...) scale_fill_manual(values = GROUP_COLORS, breaks = GROUP_ORDER, labels = GROUP_LABELS, drop = FALSE, ...)
group_shape_scale <- function(...) scale_shape_manual(values = GROUP_SHAPES, breaks = GROUP_ORDER, labels = GROUP_LABELS, drop = FALSE, ...)

# Select evenly separated decade ticks in original count units for each facet.
# Omitting 1 when zero and higher decades are visible avoids a cramped 0/1 pair.
count_axis_breaks <- function(limits) {
  limits <- sort(limits[is.finite(limits)])
  if (length(limits) < 2L || max(limits) < 0) return(numeric())
  lower <- max(0, min(limits))
  upper <- max(limits)
  if (upper == 0) return(0)
  first_power <- if (lower == 0 && upper >= 10) 1L else 0L
  last_power <- max(first_power, floor(log10(upper)))
  ticks <- c(0, 10^seq.int(first_power, last_power))
  ticks <- ticks[ticks >= lower & ticks <= upper]
  if (length(ticks) < 2L) {
    ticks <- pretty(c(lower, upper), n = 3)
    ticks <- ticks[ticks >= lower & ticks <= upper & ticks >= 0]
  }
  if (!length(ticks)) ticks <- signif(mean(c(lower, upper)), 2L)
  unique(ticks)
}

count_axis_labels <- function(values) {
  vapply(values, function(value) {
    if (!is.finite(value)) return(NA_character_)
    # Inverting log1p may return 999999.999999999 for a 1M decade tick.
    value <- signif(value, digits = 10L)
    divisor <- if (abs(value) >= 1e6) 1e6 else if (abs(value) >= 1e3) 1e3 else 1
    suffix <- if (divisor == 1e6) "M" else if (divisor == 1e3) "K" else ""
    paste0(format(signif(value / divisor, 3L), scientific = FALSE, trim = TRUE, big.mark = ""), suffix)
  }, character(1), USE.NAMES = FALSE)
}

empty_panel <- function(title, message) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No data", family = FIGURE_FONT, size = 3.5, colour = MID) +
    xlim(0, 1) + ylim(0, 1) + theme_void(base_family = FIGURE_FONT) + theme_sd()
}

write_source <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}


PANEL_EXPORT_SPEC <- data.frame(
  panel_id = c(paste0("Fig2", letters[1:4]), "Fig3", paste0("Fig4", letters[1:4])),
  folder = c(rep("Fig2_sequencing_quality", 4), "Fig3_umi_support", rep("Fig4_repertoire_overview", 4)),
  stem = c("Fig2a_read_position_quality", "Fig2b_paired_read_count_agreement", "Fig2c_per_base_n_content", "Fig2d_trimmed_read_length_distribution", "Fig3_umi_output_diversity", "Fig4a_assembled_clonotype_counts", "Fig4b_cdr3_length_summary", "Fig4c_blood_native_clonotype_counts", "Fig4d_rank_abundance"),
  width_mm = c(91.5,91.5,91.5,91.5,183,183,183,183,183),
  height_mm = c(75,75,75,75,135,78.5,24.5,30,31), stringsAsFactors=FALSE
)
yw_lw <- function(pt) pt / (0.75 * ggplot2::.pt)

export_panel <- function(plot, folder, stem, legend) {
  spec <- PANEL_EXPORT_SPEC[PANEL_EXPORT_SPEC$stem == stem, , drop = FALSE]
  if (nrow(spec) != 1L) stop("Unknown standalone panel: ", stem, call. = FALSE)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  if (!startsWith(stem, "Fig2")) {
    plot <- plot + labs(title=NULL,subtitle=NULL,caption=NULL,tag=NULL) +
      theme(text=element_text(family=FIGURE_FONT,size=6,face="plain"),
        axis.text=element_text(family=FIGURE_FONT,size=6,colour="#263238"),
        axis.text.x=element_text(family=FIGURE_FONT,size=6),
        axis.text.y=element_text(family=FIGURE_FONT,size=6),
        axis.title=element_text(family=FIGURE_FONT,size=7,face="plain"),
        axis.title.x=element_text(family=FIGURE_FONT,size=7,margin=margin(t=5)),
        axis.title.y=element_text(family=FIGURE_FONT,size=7,margin=margin(r=6)),
        axis.line=element_line(colour="black",linewidth=yw_lw(0.5)),
        axis.ticks=element_line(colour="black",linewidth=yw_lw(0.5)),
        axis.ticks.length=grid::unit(2.5,"pt"),
        panel.grid=element_blank(),
        strip.text=element_text(family=FIGURE_FONT,size=7,face="plain",margin=margin(2,1,2,1)),
        legend.title=element_text(family=FIGURE_FONT,size=6,face="plain"),
        legend.text=element_text(family=FIGURE_FONT,size=6,face="plain"),
        plot.tag=element_blank())
  } else {
    plot <- plot + labs(tag=NULL) + theme(plot.tag=element_blank())
  }
  chain_strip <- spec$panel_id %in% c("Fig4b", "Fig4c", "Fig4d")
  if (chain_strip) {
    plot <- plot + theme(
      strip.text = element_text(size = 6.5, margin = margin(0.5, 0, 0.5, 0)),
      axis.text = element_text(size = 6), axis.ticks.length = grid::unit(1.5, "pt"),
      axis.title = element_text(size = 6.5),
      axis.title.x = element_text(margin = margin(t = 1)),
      axis.title.y = element_text(margin = margin(r = 1)),
      panel.spacing.x = grid::unit(1, "mm"),
      legend.title = element_blank(), legend.position = "none",
      plot.margin = margin(1, 2, 1, 1))
    if (spec$panel_id == "Fig4b") {
      plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
  }
  width_in <- spec$width_mm / 25.4
  height_in <- spec$height_mm / 25.4
  draw <- function(open_device) {
    open_device()
    on.exit(grDevices::dev.off(), add = TRUE)
    if (spec$panel_id == "Fig4a") {
      grid::grid.newpage()
      panel <- grid::viewport(y = grid::unit(0, "mm"), width = grid::unit(177, "mm"),
                              height = grid::unit(69, "mm"),
                              just = c("centre", "bottom"))
      print(plot, newpage = FALSE, vp = panel)
      text_style <- grid::gpar(fontfamily = FIGURE_FONT, fontsize = 6, col = INK)
      label_width <- vapply(unname(GROUP_LABELS), function(label) {
        grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = text_style)), "mm", valueOnly = TRUE)
      }, numeric(1))
      entry_width <- 4.7 + label_width
      gap <- 4
      left <- grid::unit(0.5, "npc") - grid::unit((sum(entry_width) + gap) / 2, "mm")
      key_x <- left + grid::unit(c(2, entry_width[1] + gap + 2), "mm")
      y <- grid::unit(1, "npc") - grid::unit(4, "mm")
      for (i in seq_along(GROUP_ORDER)) {
        x <- key_x[i]
        colour <- unname(GROUP_COLORS[i])
        grid::grid.segments(x-grid::unit(2,"mm"), y, x+grid::unit(2,"mm"), y,
                            gp = grid::gpar(col = colour, lwd = 0.75, lty = c(1,2)[i]))
        grid::grid.points(x, y, pch = unname(GROUP_SHAPES[i]), size = grid::unit(1.4,"mm"),
                          gp = grid::gpar(col = colour, fill = colour, lwd = 0.5))
        grid::grid.text(unname(GROUP_LABELS[i]), x+grid::unit(2.7,"mm"), y, just = "left", gp = text_style)
      }
    } else if (chain_strip) {
      grid::grid.newpage()
      layout <- ggplot2::ggplotGrob(plot)
      left_axes <- layout$layout$l[grepl("^axis-l-", layout$layout$name)]
      layout$widths[left_axes] <- grid::unit(6.4, "mm")
      layout$widths[layout$layout$l[layout$layout$name == "ylab-l"]] <- grid::unit(5.6, "mm")
      panels <- sort(layout$layout$l[grepl("^panel-", layout$layout$name)])
      if (length(panels) == 7L) {
        layout <- gtable::gtable_add_cols(layout, grid::unit(2, "mm"), pos = panels[3])
      }
      grid::grid.draw(layout)
    } else print(plot)
  }
  for (format in FORMATS) {
    path <- file.path(folder, paste0(stem, ".", format))
    if (format == "pdf") {
      draw(function() grDevices::cairo_pdf(filename = path, width = width_in, height = height_in,
                                           family = FIGURE_FONT, bg = "white", antialias = "gray", fallback_resolution = 600))
    } else if (format == "svg") {
      draw(function() svglite::svglite(file = path, width = width_in, height = height_in,
                                      bg = "white", system_fonts = list(sans = FIGURE_FONT)))
    } else if (format == "tiff") {
      draw(function() grDevices::tiff(filename = path, width = width_in, height = height_in, units = "in", res = 600,
                                      type = "cairo", family = FIGURE_FONT, bg = "white", antialias = "gray", compression = "lzw"))
    } else {
      draw(function() grDevices::png(filename = path, width = width_in, height = height_in, units = "in", res = 600,
                                     type = "cairo", family = FIGURE_FONT, bg = "white", antialias = "gray"))
    }
  }
  writeLines(legend, file.path(folder, paste0(stem, "_legend.txt")), useBytes = TRUE)
}

median_iqr <- function(data, keys, value_name) {
  keep <- !is.na(data[[value_name]])
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  groups <- interaction(data[keys], drop = TRUE, lex.order = TRUE)
  pieces <- split(data, groups)
  rows <- lapply(pieces, function(z) {
    key <- z[1, keys, drop = FALSE]
    values <- z[[value_name]]
    cbind(key, n = length(values), median = median(values), q25 = unname(quantile(values, 0.25)), q75 = unname(quantile(values, 0.75)))
  })
  do.call(rbind, rows)
}

