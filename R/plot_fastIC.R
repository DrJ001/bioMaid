# ============================================================
#  plot_fastIC.R
#  ggplot2 visualisation for fastIC() output.
# ============================================================

utils::globalVariables(c(
  "global_op", "global_stab", "global_dev",
  "highlighted", "score1", "score2", "score_x", "score_y", "loads1", "loads2",
  "iclass", "iClassOP", "iClassRMSD", "CVE",
  "x_val", "y_val", "x_class", "y_class", "winner", "winning_CVE",
  "env_ord", "geno_ord", "mean_op", "label_var", "arrow_x", "arrow_y",
  "spec.var", "fitted1", "mean_iClassOP", "iclass_lbl",
  "env_label", "geno_label", "ic_label", "panel_x", "panel_y",
  "fill_val", "xend", "yend", "corr", "row_var", "col_var",
  "stab_val", "op_val", "panel_col"
))

#' @importFrom stats aggregate reshape as.formula
#' @importFrom rlang sym .data
#' @importFrom scales hue_pal
NULL

# ============================================================
#  Internal helpers  (.pfi_* prefix)
# ============================================================

# ---- Metadata extractor ------------------------------------------------

#' Extract metadata from a fastIC() result
#'
#' @param res Long-format data frame returned by fastIC().
#' @return Named list with metadata fields.
#' @noRd
.pfi_parse <- function(res) {

  nms <- names(res)

  sterm <- nms[1L]
  gterm <- nms[2L]

  load_cols   <- grep("^loads[0-9]+$",  nms, value = TRUE)
  score_cols  <- grep("^score[0-9]+$",  nms, value = TRUE)
  fitted_cols <- grep("^fitted[0-9]+$", nms, value = TRUE)

  k <- length(load_cols)

  list(
    sterm       = sterm,
    gterm       = gterm,
    k           = k,
    load_cols   = load_cols,
    score_cols  = score_cols,
    fitted_cols = fitted_cols,
    has_fast    = "global_op"   %in% nms,
    has_stab    = "global_stab" %in% nms,
    has_iclass  = "iclass"      %in% nms
  )
}

# ---- Environment-level data --------------------------------------------

#' Extract unique environment-level rows from fastIC() result
#'
#' @param res Long-format data frame returned by fastIC().
#' @param p   Metadata list from .pfi_parse().
#' @return Data frame with one row per environment.
#' @noRd
.pfi_env_data <- function(res, p) {

  keep <- c(p$sterm, p$load_cols, "spec.var",
            if (p$has_iclass) "iclass")
  keep <- intersect(keep, names(res))
  ed   <- unique(res[, keep, drop = FALSE])
  rownames(ed) <- NULL
  ed
}

# ---- Genotype-level data -----------------------------------------------

#' Extract unique genotype-level rows from fastIC() result
#'
#' Also pivots iClassOP / iClassRMSD wide when iclass is present.
#'
#' @param res Long-format data frame returned by fastIC().
#' @param p   Metadata list from .pfi_parse().
#' @return Data frame with one row per genotype.
#' @noRd
.pfi_geno_data <- function(res, p) {

  keep <- c(p$gterm, p$score_cols,
            if (p$has_fast)  "global_op",
            if (p$has_stab)  "global_stab")
  keep <- intersect(keep, names(res))
  gd   <- unique(res[, keep, drop = FALSE])
  rownames(gd) <- NULL

  # Pivot iClassOP / iClassRMSD wide: one pair of columns per iClass level
  if (p$has_iclass &&
      all(c("iclass", "iClassOP", "iClassRMSD") %in% names(res))) {

    ic_sub <- unique(res[, c(p$gterm, "iclass", "iClassOP", "iClassRMSD"),
                          drop = FALSE])
    ic_sub <- ic_sub[!is.na(ic_sub$iclass), , drop = FALSE]

    if (nrow(ic_sub) > 0L) {
      ic_levels <- levels(ic_sub$iclass)
      if (is.null(ic_levels))
        ic_levels <- sort(unique(as.character(ic_sub$iclass)))

      for (ic in ic_levels) {
        ic_rows <- ic_sub[as.character(ic_sub$iclass) == ic, , drop = FALSE]
        op_col  <- paste0("iClassOP_",   ic)
        rm_col  <- paste0("iClassRMSD_", ic)
        tmp_op  <- setNames(ic_rows[, c(p$gterm, "iClassOP")],
                            c(p$gterm, op_col))
        tmp_rm  <- setNames(ic_rows[, c(p$gterm, "iClassRMSD")],
                            c(p$gterm, rm_col))
        gd <- merge(gd, tmp_op, by = p$gterm, all.x = TRUE)
        gd <- merge(gd, tmp_rm, by = p$gterm, all.x = TRUE)
      }
    }
  }

  rownames(gd) <- NULL
  gd
}

# ---- Highlight resolver ------------------------------------------------

#' Resolve which genotypes to highlight
#'
#' @param gd          Genotype-level data frame from .pfi_geno_data().
#' @param gterm       Name of the genotype column.
#' @param highlight   User-supplied highlight argument.
#' @param n_highlight Integer; how many genotypes to highlight by default.
#' @param type        Character; the plot type being produced.
#' @return Character vector of genotype names.
#' @noRd
.pfi_highlights <- function(gd, gterm, highlight, n_highlight, type) {

  n_highlight <- max(1L, as.integer(n_highlight))
  all_genos   <- as.character(gd[[gterm]])

  # Explicit NULL — no highlights
  if (is.null(highlight)) return(character(0L))

  # User-supplied names
  if (!identical(highlight, "default") && is.character(highlight))
    return(intersect(highlight, all_genos))

  # "default" logic: branch by plot type
  if (type %in% c("iclass", "OP.pairs", "OP.variety")) {

    # Try to use mean iClassOP across classes
    ic_op_cols <- grep("^iClassOP_", names(gd), value = TRUE)

    if (length(ic_op_cols) >= 1L) {
      mat <- gd[, ic_op_cols, drop = FALSE]
      gd[["mean_iClassOP"]] <- rowMeans(mat, na.rm = TRUE)
      ord <- order(gd[["mean_iClassOP"]], decreasing = TRUE, na.last = TRUE)
      return(as.character(gd[[gterm]][head(ord, n_highlight)]))
    }

    # Fallback to global_op
    if ("global_op" %in% names(gd)) {
      ord <- order(gd[["global_op"]], decreasing = TRUE, na.last = TRUE)
      return(as.character(gd[[gterm]][head(ord, n_highlight)]))
    }

    return(as.character(head(all_genos, n_highlight)))
  }

  # "fast", "biplot", "CVE", "dev", and all others
  if ("global_op" %in% names(gd)) {

    top_op <- head(
      as.character(gd[[gterm]][
        order(gd[["global_op"]], decreasing = TRUE, na.last = TRUE)
      ]),
      n_highlight
    )

    if ("global_stab" %in% names(gd)) {
      top_stab <- head(
        as.character(gd[[gterm]][
          order(gd[["global_stab"]], decreasing = TRUE, na.last = TRUE)
        ]),
        n_highlight
      )
      return(unique(c(top_op, top_stab)))
    }

    return(top_op)
  }

  as.character(head(all_genos, n_highlight))
}

# ---- Text-repel helper -------------------------------------------------

#' Produce a text annotation layer (ggrepel if available, else geom_text)
#'
#' @param data      Data frame for the layer.
#' @param mapping   aes() mapping to use.
#' @param ...       Additional arguments forwarded to the geom.
#' @return A ggplot2 layer.
#' @noRd
.pfi_text_layer <- function(data, mapping, ...) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = data, mapping = mapping,
                             size = 2.8, max.overlaps = 20L,
                             segment.colour = "grey60",
                             segment.size   = 0.3,
                             show.legend    = FALSE,
                             ...)
  } else {
    ggplot2::geom_text(data = data, mapping = mapping,
                       size = 2.8, vjust = -0.6,
                       show.legend = FALSE,
                       ...)
  }
}

# ---- Qualitative colour palette ----------------------------------------

#' Return n qualitative colours
#' @noRd
.pfi_pal <- function(n) {
  if (n <= 0L) return(character(0L))
  scales::hue_pal()(n)
}

# ---- Heatmap builder (shared by "cve" and "dev") -----------------------

#' Build a GEI heatmap (CVE or dev)
#'
#' @param res       Long-format fastIC() result.
#' @param p         Metadata list from .pfi_parse().
#' @param fill_col  Character; "CVE" or "dev".
#' @param fill_lab  Character; legend title.
#' @param title     Character; plot title.
#' @param theme     ggplot2 theme object.
#' @return ggplot object.
#' @noRd
.pfi_heatmap <- function(res, p, fill_col, fill_lab, title, theme) {

  df <- res[, c(p$sterm, p$gterm, fill_col,
                if (p$has_iclass) "iclass",
                if ("loads1" %in% names(res)) "loads1"),
            drop = FALSE]

  # ---- Order environments
  env_levs <- if (p$has_iclass && "iclass" %in% names(df)) {
    df2 <- unique(df[, c(p$sterm,
                          if ("loads1" %in% names(df)) "loads1",
                          "iclass"), drop = FALSE])
    df2 <- df2[order(df2$iclass,
                     if ("loads1" %in% names(df2)) -df2$loads1
                     else seq_len(nrow(df2))), ]
    as.character(df2[[p$sterm]])
  } else if ("loads1" %in% names(df)) {
    df2 <- unique(df[, c(p$sterm, "loads1"), drop = FALSE])
    df2 <- df2[order(-df2$loads1), ]
    as.character(df2[[p$sterm]])
  } else {
    sort(unique(as.character(df[[p$sterm]])))
  }

  # ---- Order genotypes
  gd       <- .pfi_geno_data(res, p)
  geno_ord <- if ("global_op" %in% names(gd)) {
    as.character(gd[[p$gterm]][order(gd$global_op, decreasing = TRUE, na.last = TRUE)])
  } else {
    sort(unique(as.character(df[[p$gterm]])))
  }

  df[[p$sterm]] <- factor(df[[p$sterm]], levels = env_levs)
  df[[p$gterm]] <- factor(df[[p$gterm]], levels = rev(geno_ord))

  names(df)[names(df) == fill_col] <- "fill_val"

  # ---- Build plot
  if (p$has_iclass && "iclass" %in% names(df)) {

    p_out <- ggplot2::ggplot(df,
               ggplot2::aes(x = !!rlang::sym(p$sterm),
                            y = !!rlang::sym(p$gterm),
                            fill = fill_val)) +
      ggplot2::geom_tile(colour = NA) +
      ggplot2::scale_fill_gradient2(
        low      = "#B2182B",
        mid      = "white",
        high     = "#2166AC",
        midpoint = 0,
        name     = fill_lab
      ) +
      ggplot2::facet_grid(
        stats::as.formula(paste(". ~", "iclass")),
        scales = "free_x",
        space  = "free_x"
      ) +
      ggplot2::labs(
        x     = p$sterm,
        y     = p$gterm,
        title = title
      ) +
      theme +
      ggplot2::theme(
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                                  size = 7),
        axis.text.y      = ggplot2::element_text(size = 7),
        strip.background = ggplot2::element_rect(fill = "grey88",
                                                  colour = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold"),
        panel.spacing    = ggplot2::unit(0.3, "lines")
      )

  } else {

    p_out <- ggplot2::ggplot(df,
               ggplot2::aes(x = !!rlang::sym(p$sterm),
                            y = !!rlang::sym(p$gterm),
                            fill = fill_val)) +
      ggplot2::geom_tile(colour = NA) +
      ggplot2::scale_fill_gradient2(
        low      = "#B2182B",
        mid      = "white",
        high     = "#2166AC",
        midpoint = 0,
        name     = fill_lab
      ) +
      ggplot2::labs(
        x     = p$sterm,
        y     = p$gterm,
        title = title
      ) +
      theme +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = ggplot2::element_text(size = 7)
      )
  }

  # ---- Wide data for return_data
  wide_df <- stats::reshape(
    df[, c(p$sterm, p$gterm, "fill_val"), drop = FALSE],
    idvar     = p$gterm,
    timevar   = p$sterm,
    direction = "wide",
    v.names   = "fill_val"
  )
  names(wide_df) <- sub("^fill_val\\.", "", names(wide_df))

  list(plot = p_out, data = wide_df)
}


# ============================================================
#  Plot-type builders  (.pfi_plot_* prefix)
# ============================================================

# ---- "fast" ------------------------------------------------------------

#' @noRd
.pfi_plot_fast <- function(res, p, hl_names, theme) {

  gd <- .pfi_geno_data(res, p)

  gd$highlighted <- gd[[p$gterm]] %in% hl_names
  n_hl           <- length(hl_names)
  hl_cols        <- if (n_hl > 0L) .pfi_pal(n_hl) else character(0L)

  gd_hl  <- gd[gd$highlighted,  , drop = FALSE]
  gd_reg <- gd[!gd$highlighted, , drop = FALSE]

  mean_op   <- mean(gd$global_op,   na.rm = TRUE)
  mean_stab <- mean(gd$global_stab, na.rm = TRUE)

  x_rng <- range(gd$global_op,   na.rm = TRUE)
  y_rng <- range(gd$global_stab, na.rm = TRUE)

  x_lo <- x_rng[1L]; x_hi <- x_rng[2L]
  y_lo <- y_rng[1L]; y_hi <- y_rng[2L]

  # Small inset padding
  px <- 0.03 * diff(x_rng)
  py <- 0.03 * diff(y_rng)

  plt <- ggplot2::ggplot(gd, ggplot2::aes(x = global_op, y = global_stab)) +
    ggplot2::geom_hline(yintercept = mean_stab,
                        linetype = "dashed", colour = "grey60") +
    ggplot2::geom_vline(xintercept = mean_op,
                        linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(data = gd_reg, size = 2, colour = "grey60") +
    ggplot2::annotate("text",
      x      = x_hi - px, y = y_lo + py,
      label  = "Broadly\nadapted",
      hjust  = 1, vjust = 0, alpha = 0.4, size = 3.5
    ) +
    ggplot2::annotate("text",
      x      = x_hi - px, y = y_hi - py,
      label  = "Responsive",
      hjust  = 1, vjust = 1, alpha = 0.4, size = 3.5
    ) +
    ggplot2::annotate("text",
      x      = x_lo + px, y = y_lo + py,
      label  = "Poor\n& stable",
      hjust  = 0, vjust = 0, alpha = 0.4, size = 3.5
    ) +
    ggplot2::annotate("text",
      x      = x_lo + px, y = y_hi - py,
      label  = "Poor\n& unstable",
      hjust  = 0, vjust = 1, alpha = 0.4, size = 3.5
    )

  if (n_hl > 0L) {
    gd_hl[[p$gterm]] <- factor(gd_hl[[p$gterm]], levels = hl_names)
    plt <- plt +
      ggplot2::geom_point(
        data    = gd_hl,
        mapping = ggplot2::aes(colour = !!rlang::sym(p$gterm)),
        size    = 2.5
      ) +
      ggplot2::scale_colour_manual(values = hl_cols, name = p$gterm) +
      .pfi_text_layer(
        data    = gd_hl,
        mapping = ggplot2::aes(label = !!rlang::sym(p$gterm),
                               colour = !!rlang::sym(p$gterm))
      )
  }

  plt <- plt +
    ggplot2::labs(
      x     = "Global Overall Performance (global_op)",
      y     = "Global Stability (global_stab)",
      title = "FAST: Global OP vs Stability"
    ) +
    theme

  gd$highlighted <- gd[[p$gterm]] %in% hl_names
  list(plot = plt, data = gd)
}

# ---- "biplot" ----------------------------------------------------------

#' @noRd
.pfi_plot_biplot <- function(res, p, hl_names, theme,
                             biplot_factors = c(1L, 2L)) {

  gd <- .pfi_geno_data(res, p)
  ed <- .pfi_env_data(res, p)

  gd$highlighted <- gd[[p$gterm]] %in% hl_names
  n_hl           <- length(hl_names)

  # Select the requested factor axes
  fa <- biplot_factors
  sc_x <- paste0("score", fa[1L]);  sc_y <- paste0("score", fa[2L])
  ld_x <- paste0("loads", fa[1L]);  ld_y <- paste0("loads", fa[2L])

  # Scaling: fit arrows into score space
  score_range <- max(abs(c(gd[[sc_x]], gd[[sc_y]])), na.rm = TRUE)
  load_range  <- max(abs(c(ed[[ld_x]], ed[[ld_y]])), na.rm = TRUE)

  arrow_scale  <- if (load_range > 0) 0.85 * score_range / load_range else 1
  ed_sc        <- ed
  ed_sc$loads1 <- ed[[ld_x]] * arrow_scale
  ed_sc$loads2 <- ed[[ld_y]] * arrow_scale

  # Alias the chosen score axes so downstream aes() can use fixed names
  gd$score_x <- gd[[sc_x]]
  gd$score_y <- gd[[sc_y]]

  gd_hl  <- gd[gd$highlighted,  , drop = FALSE]
  gd_reg <- gd[!gd$highlighted, , drop = FALSE]
  hl_cols <- if (n_hl > 0L) .pfi_pal(n_hl) else character(0L)

  # Arrow colour: iclass (if present) or grey40
  if (p$has_iclass && "iclass" %in% names(ed_sc)) {
    ic_levels <- levels(ed_sc$iclass)
    if (is.null(ic_levels))
      ic_levels <- sort(unique(as.character(ed_sc$iclass)))
    ic_cols <- .pfi_pal(length(ic_levels))
    names(ic_cols) <- ic_levels
    arrow_col_arg <- TRUE
  } else {
    arrow_col_arg <- FALSE
  }

  plt <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80")

  if (arrow_col_arg) {
    plt <- plt +
      ggplot2::geom_segment(
        data    = ed_sc,
        mapping = ggplot2::aes(
          x      = 0, y = 0,
          xend   = loads1, yend = loads2,
          colour = iclass
        ),
        arrow     = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
        linewidth = 0.7
      ) +
      ggplot2::scale_colour_manual(
        values = ic_cols,
        name   = "iClass",
        guide  = ggplot2::guide_legend(order = 1L)
      ) +
      .pfi_text_layer(
        data    = ed_sc,
        mapping = ggplot2::aes(x = loads1, y = loads2,
                               label = !!rlang::sym(p$sterm),
                               colour = iclass)
      )
  } else {
    plt <- plt +
      ggplot2::geom_segment(
        data    = ed_sc,
        mapping = ggplot2::aes(x = 0, y = 0, xend = loads1, yend = loads2),
        colour   = "grey40",
        arrow    = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
        linewidth = 0.7
      ) +
      .pfi_text_layer(
        data    = ed_sc,
        mapping = ggplot2::aes(x = loads1, y = loads2,
                               label = !!rlang::sym(p$sterm)),
        colour  = "grey30"
      )
  }

  plt <- plt +
    ggplot2::geom_point(
      data    = gd_reg,
      mapping = ggplot2::aes(x = score_x, y = score_y),
      colour  = "grey70",
      size    = 2
    )

  if (n_hl > 0L) {
    gd_hl[[p$gterm]] <- factor(gd_hl[[p$gterm]], levels = hl_names)
    plt <- plt +
      ggplot2::geom_point(
        data    = gd_hl,
        mapping = ggplot2::aes(x = score_x, y = score_y,
                               fill = !!rlang::sym(p$gterm)),
        shape   = 21,
        colour  = "white",
        size    = 3
      ) +
      ggplot2::scale_fill_manual(
        values = hl_cols,
        name   = p$gterm,
        guide  = ggplot2::guide_legend(order = 2L)
      ) +
      .pfi_text_layer(
        data    = gd_hl,
        mapping = ggplot2::aes(x = score_x, y = score_y,
                               label = !!rlang::sym(p$gterm))
      )
  }

  plt <- plt +
    ggplot2::labs(
      x     = paste0("Factor ", fa[1L], " score / loading"),
      y     = paste0("Factor ", fa[2L], " score / loading"),
      title = paste0("FA biplot: Genotypes (points) and Environments (arrows)",
                     if (!identical(fa, c(1L, 2L)))
                       paste0(" \u2014 Factors ", fa[1L], " & ", fa[2L])
                     else "")
    ) +
    theme

  list(plot = plt, data = list(geno = gd, envs = ed_sc))
}

# ---- "CVE" -------------------------------------------------------------

#' @noRd
.pfi_plot_CVE <- function(res, p, theme) {
  .pfi_heatmap(
    res      = res,
    p        = p,
    fill_col = "CVE",
    fill_lab = "CVE",
    title    = "Common Variety Effect (CVE) \u2014 Genotype x Environment",
    theme    = theme
  )
}

# ---- "iclass" ----------------------------------------------------------

#' @noRd
.pfi_plot_iclass <- function(res, p, hl_names, theme) {

  ic_df <- unique(res[, c(p$gterm, "iclass", "iClassOP", "iClassRMSD"),
                       drop = FALSE])
  ic_df <- ic_df[!is.na(ic_df$iclass), , drop = FALSE]

  ic_df$highlighted <- ic_df[[p$gterm]] %in% hl_names
  n_hl     <- length(hl_names)
  hl_cols  <- if (n_hl > 0L) .pfi_pal(n_hl) else character(0L)

  # Compute per-class mean iClassOP for the vertical reference line
  class_means <- stats::aggregate(
    iClassOP ~ iclass,
    data = ic_df,
    FUN  = mean, na.rm = TRUE
  )
  names(class_means)[2L] <- "mean_ic_op"

  ic_df <- merge(ic_df, class_means, by = "iclass", all.x = TRUE)

  # Count environments per iClass for facet labeller
  env_df  <- .pfi_env_data(res, p)
  if ("iclass" %in% names(env_df)) {
    ic_n   <- table(as.character(env_df$iclass))
    ic_lab <- function(x) {
      paste0(x, " (n=", ic_n[x], " envs)")
    }
  } else {
    ic_lab <- function(x) x
  }
  lbl_fn <- ggplot2::as_labeller(ic_lab)

  ic_reg <- ic_df[!ic_df$highlighted, , drop = FALSE]
  ic_hl  <- ic_df[ic_df$highlighted,  , drop = FALSE]

  plt <- ggplot2::ggplot(ic_df,
    ggplot2::aes(x = iClassOP, y = iClassRMSD)) +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = "dashed", colour = "grey60") +
    ggplot2::geom_vline(
      data    = unique(ic_df[, c("iclass", "mean_ic_op"), drop = FALSE]),
      mapping = ggplot2::aes(xintercept = mean_ic_op),
      linetype = "dashed", colour = "grey60"
    ) +
    ggplot2::geom_point(data = ic_reg, colour = "grey70", size = 2)

  if (n_hl > 0L) {
    ic_hl[[p$gterm]] <- factor(ic_hl[[p$gterm]], levels = hl_names)
    plt <- plt +
      ggplot2::geom_point(
        data    = ic_hl,
        mapping = ggplot2::aes(colour = !!rlang::sym(p$gterm)),
        size    = 2.5
      ) +
      ggplot2::scale_colour_manual(values = hl_cols, name = p$gterm) +
      .pfi_text_layer(
        data    = ic_hl,
        mapping = ggplot2::aes(label = !!rlang::sym(p$gterm),
                               colour = !!rlang::sym(p$gterm))
      )
  }

  plt <- plt +
    ggplot2::facet_wrap(
      stats::as.formula("~ iclass"),
      scales   = "free",
      labeller = lbl_fn
    ) +
    ggplot2::labs(
      x     = "iClassOP",
      y     = "iClassRMSD",
      title = "iClass: Within-class OP vs RMSD"
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey88",
                                                colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold")
    )

  list(plot = plt, data = ic_df)
}

# ---- "OP.variety" -------------------------------------------------------

#' @noRd
.pfi_plot_op_variety <- function(res, p, hl_names, theme) {

  traj_df <- unique(res[, c(p$gterm, "iclass", "iClassOP"), drop = FALSE])
  traj_df <- traj_df[!is.na(traj_df$iclass), , drop = FALSE]

  # Preserve iclass level order
  if (!is.factor(traj_df$iclass)) {
    if (!is.null(levels(res$iclass))) {
      traj_df$iclass <- factor(traj_df$iclass, levels = levels(res$iclass))
    } else {
      traj_df$iclass <- factor(traj_df$iclass,
                                levels = sort(unique(traj_df$iclass)))
    }
  }

  traj_df$highlighted <- traj_df[[p$gterm]] %in% hl_names
  n_hl    <- length(hl_names)
  hl_cols <- if (n_hl > 0L) .pfi_pal(n_hl) else character(0L)

  bg_df <- traj_df[!traj_df$highlighted, , drop = FALSE]
  hl_df <- traj_df[traj_df$highlighted,  , drop = FALSE]

  # Rightmost iClass level for label placement
  ic_levels <- levels(traj_df$iclass)
  last_ic   <- ic_levels[length(ic_levels)]
  label_df  <- hl_df[as.character(hl_df$iclass) == last_ic, , drop = FALSE]

  plt <- ggplot2::ggplot(
    traj_df,
    ggplot2::aes(x = iclass, y = iClassOP,
                 group = !!rlang::sym(p$gterm))
  ) +
    ggplot2::geom_line(
      data      = bg_df,
      colour    = "grey85",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      data   = bg_df,
      colour = "grey85",
      size   = 1
    )

  if (n_hl > 0L) {
    hl_df[[p$gterm]] <- factor(hl_df[[p$gterm]], levels = hl_names)
    if (nrow(label_df) > 0L)
      label_df[[p$gterm]] <- factor(label_df[[p$gterm]], levels = hl_names)

    plt <- plt +
      ggplot2::geom_line(
        data      = hl_df,
        mapping   = ggplot2::aes(colour = !!rlang::sym(p$gterm)),
        linewidth = 1.1
      ) +
      ggplot2::geom_point(
        data    = hl_df,
        mapping = ggplot2::aes(colour = !!rlang::sym(p$gterm)),
        size    = 3
      ) +
      ggplot2::scale_colour_manual(values = hl_cols, name = p$gterm)

    if (nrow(label_df) > 0L) {
      plt <- plt +
        .pfi_text_layer(
          data    = label_df,
          mapping = ggplot2::aes(label  = !!rlang::sym(p$gterm),
                                 colour = !!rlang::sym(p$gterm)),
          nudge_x = 0.2
        )
    }
  }

  plt <- plt +
    ggplot2::labs(
      x     = "iClass",
      y     = "iClassOP",
      title = "iClassOP across interaction classes (variety trajectories)"
    ) +
    theme

  list(plot = plt, data = traj_df)
}

# ---- "OP.pairs" --------------------------------------------------------

#' @noRd
.pfi_plot_op_pairs <- function(res, p, hl_names, theme) {

  traj_df <- unique(res[, c(p$gterm, "iclass", "iClassOP"), drop = FALSE])
  traj_df <- traj_df[!is.na(traj_df$iclass), , drop = FALSE]

  ic_levels <- if (!is.null(levels(res$iclass))) {
    levels(res$iclass)
  } else {
    sort(unique(as.character(traj_df$iclass)))
  }

  n_ic <- length(ic_levels)

  if (n_ic < 2L) {
    warning("plot_fastIC(): 'pairs' plot requires >= 2 iClass levels; ",
            "only ", n_ic, " found. Returning empty plot.")
    plt <- ggplot2::ggplot() +
      ggplot2::labs(title = "OP.pairs: insufficient iClass levels") +
      theme
    return(list(plot = plt, data = traj_df))
  }

  # Pivot wide: one iClassOP column per class
  wide <- as.data.frame(matrix(NA_real_,
    nrow = length(unique(traj_df[[p$gterm]])),
    ncol = 1L + n_ic
  ))
  genotypes <- sort(unique(as.character(traj_df[[p$gterm]])))
  wide[[1L]] <- genotypes
  names(wide) <- c(p$gterm, paste0("iClassOP_", ic_levels))

  for (ic in ic_levels) {
    sub <- traj_df[as.character(traj_df$iclass) == ic, , drop = FALSE]
    m   <- match(genotypes, as.character(sub[[p$gterm]]))
    wide[[paste0("iClassOP_", ic)]] <- sub$iClassOP[m]
  }

  wide$highlighted <- wide[[p$gterm]] %in% hl_names
  n_hl    <- length(hl_names)
  hl_cols <- if (n_hl > 0L) .pfi_pal(n_hl) else character(0L)

  # Build lower-triangular panel data frames
  pairs_rows <- vector("list", n_ic * (n_ic - 1L) / 2L)
  pair_idx   <- 1L

  for (i in seq(2L, n_ic)) {
    for (j in seq(1L, i - 1L)) {
      col_x  <- paste0("iClassOP_", ic_levels[j])
      col_y  <- paste0("iClassOP_", ic_levels[i])
      tmp    <- data.frame(
        x_class     = ic_levels[j],
        y_class     = ic_levels[i],
        x_val       = wide[[col_x]],
        y_val       = wide[[col_y]],
        highlighted = wide$highlighted,
        stringsAsFactors = FALSE
      )
      tmp[[p$gterm]] <- wide[[p$gterm]]
      pairs_rows[[pair_idx]] <- tmp
      pair_idx <- pair_idx + 1L
    }
  }

  pairs_df              <- do.call(rbind, pairs_rows)
  pairs_df$x_class      <- factor(pairs_df$x_class, levels = ic_levels)
  pairs_df$y_class      <- factor(pairs_df$y_class, levels = ic_levels)

  pairs_reg <- pairs_df[!pairs_df$highlighted, , drop = FALSE]
  pairs_hl  <- pairs_df[pairs_df$highlighted,  , drop = FALSE]

  # Try patchwork
  use_pw <- requireNamespace("patchwork", quietly = TRUE)

  if (use_pw) {
    # Build individual panel plots and assemble.
    # Layout is (n_ic-1) x (n_ic-1); lower-triangular panels hold real data,
    # diagonal and upper-triangular are plot_spacer() (true empty space).
    # y-axis shown only in leftmost column (j_idx == 1).
    # x-axis shown only in bottom row    (i_idx == n_ic).
    panel_list <- vector("list", (n_ic - 1L) * (n_ic - 1L))

    for (i_idx in seq(2L, n_ic)) {
      for (j_idx in seq(1L, n_ic - 1L)) {

        p_idx <- (i_idx - 2L) * (n_ic - 1L) + j_idx

        if (j_idx < i_idx) {
          # ---- Lower-triangular: real data panel --------------------------
          sub_reg <- pairs_reg[
            pairs_reg$y_class == ic_levels[i_idx] &
            pairs_reg$x_class == ic_levels[j_idx], , drop = FALSE
          ]
          sub_hl <- pairs_hl[
            pairs_hl$y_class == ic_levels[i_idx] &
            pairs_hl$x_class == ic_levels[j_idx], , drop = FALSE
          ]

          pp <- ggplot2::ggplot() +
            ggplot2::geom_abline(
              slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50"
            ) +
            ggplot2::geom_point(
              data    = sub_reg,
              mapping = ggplot2::aes(x = x_val, y = y_val),
              colour  = "grey70",
              size    = 1.5
            )

          if (n_hl > 0L && nrow(sub_hl) > 0L) {
            sub_hl[[p$gterm]] <- factor(sub_hl[[p$gterm]], levels = hl_names)
            pp <- pp +
              ggplot2::geom_point(
                data    = sub_hl,
                mapping = ggplot2::aes(x = x_val, y = y_val,
                                       colour = !!rlang::sym(p$gterm)),
                size    = 2.5
              ) +
              ggplot2::scale_colour_manual(values = hl_cols, name = p$gterm) +
              .pfi_text_layer(
                data    = sub_hl,
                mapping = ggplot2::aes(x = x_val, y = y_val,
                                       label = !!rlang::sym(p$gterm),
                                       colour = !!rlang::sym(p$gterm))
              )
          }

          # Axis labels only on border panels
          x_lab <- if (i_idx == n_ic)  ic_levels[j_idx]  else ""
          y_lab <- if (j_idx == 1L)    ic_levels[i_idx]  else ""

          pp <- pp +
            ggplot2::labs(x = x_lab, y = y_lab) +
            theme

          # Strip axis text/ticks from non-border panels to maximise space
          if (j_idx > 1L)
            pp <- pp + ggplot2::theme(
              axis.text.y  = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank()
            )
          if (i_idx < n_ic)
            pp <- pp + ggplot2::theme(
              axis.text.x  = ggplot2::element_blank(),
              axis.ticks.x = ggplot2::element_blank()
            )

          panel_list[[p_idx]] <- pp

        } else {
          # ---- Diagonal or upper-triangular: true empty space -------------
          panel_list[[p_idx]] <- patchwork::plot_spacer()
        }
      }
    }

    plt <- patchwork::wrap_plots(panel_list, ncol = n_ic - 1L) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = "Pairs plot of iClassOP across interaction classes"
      )
    plt <- plt & ggplot2::theme(legend.position = "right")

  } else {
    # Fallback: facet_grid on lower-triangular subset
    message("plot_fastIC(): 'patchwork' not available; ",
            "using facet_grid for pairs plot (lower triangle only).")

    plt <- ggplot2::ggplot(pairs_df,
      ggplot2::aes(x = x_val, y = y_val)) +
      ggplot2::geom_abline(
        slope = 1, intercept = 0,
        linetype = "dashed", colour = "grey50"
      ) +
      ggplot2::geom_point(
        data    = pairs_reg,
        colour  = "grey70",
        size    = 1.5
      )

    if (n_hl > 0L && nrow(pairs_hl) > 0L) {
      pairs_hl[[p$gterm]] <- factor(pairs_hl[[p$gterm]], levels = hl_names)
      plt <- plt +
        ggplot2::geom_point(
          data    = pairs_hl,
          mapping = ggplot2::aes(colour = !!rlang::sym(p$gterm)),
          size    = 2.5
        ) +
        ggplot2::scale_colour_manual(values = hl_cols, name = p$gterm) +
        .pfi_text_layer(
          data    = pairs_hl,
          mapping = ggplot2::aes(label  = !!rlang::sym(p$gterm),
                                 colour = !!rlang::sym(p$gterm))
        )
    }

    plt <- plt +
      ggplot2::facet_grid(
        y_class ~ x_class,
        scales = "free"
      ) +
      ggplot2::labs(
        x     = "iClassOP (column class)",
        y     = "iClassOP (row class)",
        title = "Pairs plot of iClassOP across interaction classes"
      ) +
      theme +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey88",
                                                  colour = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold")
      )
  }

  list(plot = plt, data = pairs_df)
}

# ============================================================
#  Public function
# ============================================================

#' Plot Output from fastIC()
#'
#' @description
#' Produces one of several ggplot2 visualisations from a [fastIC()] result
#' object.  The function dispatches to a specialised plot builder depending
#' on the `type` argument and returns either a `ggplot` object or a list
#' containing the plot and the underlying tidy data.
#'
#' @details
#' ## Plot types
#' \describe{
#'   \item{`"fast"`}{Scatter of global Overall Performance (`global_op`) vs
#'     global stability (`global_stab`) with quadrant annotations.}
#'   \item{`"biplot"`}{Factor-analysis biplot with genotype score points and
#'     environment loading arrows.  Default axes are Factors 1 \& 2; use
#'     `biplot_factors` to display any two factor axes.  Requires k
#'     \eqn{\geq} 2.}
#'   \item{`"CVE"`}{Diverging-colour heatmap of the Common Variety Effect
#'     (genotype \eqn{\times} environment), ordered by iClass then
#'     first-factor loading.}
#'   \item{`"iclass"`}{Scatter of within-class OP vs within-class RMSD,
#'     faceted by iClass.}
#'   \item{`"OP.pairs"`}{Lower-triangular pairs plot of iClassOP values across
#'     all iClass levels.  Uses \pkg{patchwork} if available; falls back to
#'     `facet_grid`.  Requires \eqn{\geq} 2 iClass levels.}
#'   \item{`"OP.variety"`}{Line plot of iClassOP across ordered interaction
#'     classes for highlighted varieties, with all others as a grey background.}
#' }
#'
#' ## Highlighting
#' When `highlight = "default"` the function automatically selects genotypes
#' to annotate based on the plot type:
#' \itemize{
#'   \item For `"fast"`, `"biplot"`, `"CVE"`: top `n_highlight` by
#'     `global_op` \emph{and} top `n_highlight` by instability (`global_stab`),
#'     deduplicated.
#'   \item For `"iclass"`, `"OP.pairs"`, `"OP.variety"`: top `n_highlight` by
#'     mean iClassOP across classes (or by `global_op` as fallback).
#' }
#' Supply `highlight = NULL` to suppress all annotation, or a character vector
#' of variety names to annotate specific genotypes.
#'
#' ## Dependencies
#' Annotation labels use \pkg{ggrepel} when available (non-overlapping), with
#' a fallback to `geom_text`.  The `"pairs"` plot uses \pkg{patchwork} when
#' available for individual panel assembly.  Colour palettes come from
#' \pkg{scales}.
#'
#' @param res            A long-format data frame returned by [fastIC()].
#' @param type           Character; the visualisation to produce. One of
#'   `"fast"` (default), `"biplot"`, `"CVE"`, `"iclass"`,
#'   `"OP.pairs"`, `"OP.variety"`.  May be abbreviated.  The Variance
#'   Accounted For chart now lives in [plot_faSummary()] with `type = "VAF"`.
#' @param highlight      Controls variety annotation. One of `"default"`
#'   (automatic; see Details), a character vector of variety names, or `NULL`
#'   (no annotation).  Default `"default"`.
#' @param n_highlight    Positive integer. Maximum number of genotypes to
#'   highlight automatically when `highlight = "default"`. Default `3L`.
#' @param biplot_factors Integer vector of length 2 giving the FA factor
#'   indices to use as the x- and y-axes of the biplot. Default `c(1L, 2L)`.
#'   Both values must be distinct and within `1:k` where `k` is the number of
#'   FA factors in `res`.  Useful when `ic.num >= 3` and iClass separation
#'   involves a third factor that is invisible in the default Factor 1 vs
#'   Factor 2 view — e.g. `biplot_factors = c(1L, 3L)` to reveal Factor 3
#'   separation.  Ignored for all plot types other than `"biplot"`.
#' @param theme          A complete ggplot2 theme object. Default
#'   [ggplot2::theme_bw()].
#' @param return_data    Logical. If `TRUE`, returns a list with elements
#'   `plot` (the ggplot object) and `data` (the underlying tidy data). Default
#'   `FALSE`.
#' @param ...         Reserved for future use; currently ignored.
#'
#' @return When `return_data = FALSE` (default): a `ggplot` object.
#'   When `return_data = TRUE`: a named list with elements `plot` and `data`.
#'   The structure of `data` is type-specific — see the Details section for
#'   each plot type.
#'
#' @seealso [fastIC()], [faSummary()], [plot_faSummary()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#' res <- fastIC(model, ic.num = 2L)
#'
#' # FAST scatter
#' plot_fastIC(res, type = "fast")
#'
#' # FA biplot with custom highlights
#' plot_fastIC(res, type = "biplot", highlight = c("Var01", "Var22"))
#'
#' # CVE heatmap, retrieve data
#' out <- plot_fastIC(res, type = "CVE", return_data = TRUE)
#' head(out$data)
#'
#' # iClassOP trajectory for top 5 varieties
#' plot_fastIC(res, type = "OP.variety", n_highlight = 5L)
#'
#' # Lower-triangular iClassOP pairs plot
#' plot_fastIC(res, type = "OP.pairs")
#' }
#'
#' @export
plot_fastIC <- function(res,
                        type           = c("fast", "biplot", "CVE",
                                           "iclass", "OP.pairs", "OP.variety"),
                        highlight      = "default",
                        n_highlight    = 3L,
                        biplot_factors = c(1L, 2L),
                        theme          = ggplot2::theme_bw(),
                        return_data    = FALSE,
                        ...) {

  # ------------------------------------------------------------------ #
  #  Validate inputs                                                     #
  # ------------------------------------------------------------------ #
  if (!is.data.frame(res) || ncol(res) < 4L)
    stop("plot_fastIC(): 'res' must be a data.frame with at least 4 columns. ",
         "Was it produced by fastIC()?")

  type        <- match.arg(type)
  n_highlight <- max(1L, as.integer(n_highlight[1L]))

  biplot_factors <- as.integer(biplot_factors)
  if (length(biplot_factors) != 2L || anyNA(biplot_factors) ||
      biplot_factors[1L] == biplot_factors[2L])
    stop("plot_fastIC(): 'biplot_factors' must be a length-2 integer vector ",
         "of distinct factor indices, e.g. c(1L, 2L).")

  if (!inherits(theme, "theme"))
    stop("plot_fastIC(): 'theme' must be a ggplot2 theme object, ",
         "e.g. ggplot2::theme_bw().")

  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("plot_fastIC(): 'return_data' must be a single logical value.")

  # Parse metadata
  p <- .pfi_parse(res)

  # ---- Type-specific precondition checks
  if (type == "fast") {
    if (!p$has_fast)
      stop("plot_fastIC(): 'global_op' column not found. ",
           "Was 'res' produced by fastIC()?")
    if (!p$has_stab)
      stop("plot_fastIC(): 'global_stab' column not found; ",
           "re-run fastIC() with k > 1.")
  }

  if (type == "biplot") {
    if (p$k < 2L)
      stop("plot_fastIC(): type = 'biplot' requires k >= 2 factors; ",
           "only ", p$k, " factor(s) found.")
    if (any(biplot_factors < 1L) || any(biplot_factors > p$k))
      stop("plot_fastIC(): 'biplot_factors' values must be between 1 and k = ",
           p$k, "; got c(", biplot_factors[1L], ", ", biplot_factors[2L], ").")
  }

  if (type %in% c("iclass", "OP.variety", "OP.pairs")) {
    if (!p$has_iclass)
      stop("plot_fastIC(): 'iclass' column not found. ",
           "Was 'res' produced by fastIC()?")
  }

  # ------------------------------------------------------------------ #
  #  Resolve highlights                                                  #
  # ------------------------------------------------------------------ #
  gd       <- .pfi_geno_data(res, p)
  hl_names <- .pfi_highlights(gd, p$gterm, highlight, n_highlight, type)

  # ------------------------------------------------------------------ #
  #  Dispatch to plot builders                                           #
  # ------------------------------------------------------------------ #
  result <- switch(type,

    fast       = .pfi_plot_fast(      res, p, hl_names, theme),
    biplot     = .pfi_plot_biplot(    res, p, hl_names, theme, biplot_factors),
    CVE        = .pfi_plot_CVE(       res, p,           theme),
    iclass     = .pfi_plot_iclass(    res, p, hl_names, theme),
    OP.pairs   = .pfi_plot_op_pairs(  res, p, hl_names, theme),
    OP.variety = .pfi_plot_op_variety(res, p, hl_names, theme)
  )

  # ------------------------------------------------------------------ #
  #  Return                                                              #
  # ------------------------------------------------------------------ #
  if (return_data) {
    result   # already a list(plot=, data=)
  } else {
    result$plot
  }
}
