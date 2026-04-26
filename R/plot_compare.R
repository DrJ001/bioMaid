# ============================================================
#  plot_compare.R
#  ggplot2 visualisations for compare() output.
#
#  Three plot types:
#   "dotplot"  — sorted predicted values with shaded criterion band
#   "letters"  — compact letter display (CLD) on a sorted dot plot
#   "heatmap"  — pairwise difference matrix coloured by significance
#
#  Design:
#   - res already carries the criterion as a named column (HSD / LSD /
#     Bonferroni), so no re-computation is needed.
#   - Group detection is automatic: columns that are constant within
#     each block of rows sharing the same criterion value are "by"
#     columns; the remaining factor column is the comparison variable.
#   - Multi-factor by groups produce colon-joined labels and are
#     treated as a single facet variable.
#   - All plot types support return_data = TRUE for bespoke ggplot2
#     customisation.
#   - reference (dotplot only): name of a variety to anchor the
#     criterion band. NULL = top-ranked (default). When set, the band
#     spans [ref - crit, ref + crit] and points are coloured three ways:
#     green (significantly better), blue (not sig.), red (worse).
#   - interactive: when TRUE, converts the ggplot to a plotly object
#     via ggplotly() with hover tooltips. Requires the plotly package.
# ============================================================


# ---- pc_interactive wrapper class -----------------------------------------
#
# When interactive = TRUE, plot_compare() returns a pc_interactive object
# rather than converting to plotly immediately.  This allows the user to
# extend the plot with ggplot2's + operator BEFORE the plotly conversion:
#
#   plot_compare(res, interactive = TRUE) + ggtitle("My title")
#
# The conversion to plotly happens lazily when the object is printed.

## Wrap a ggplot in the pc_interactive class.
#' @noRd
.pc_make_interactive <- function(p) {
  structure(list(plot = p), class = "pc_interactive")
}

## Support ggplot2's + operator on pc_interactive objects.
## Adds the layer to the internal ggplot, then re-wraps.
#' @export
`+.pc_interactive` <- function(e1, e2) {
  e1$plot <- e1$plot + e2
  e1
}

## Auto-convert to plotly on print.
#' @export
print.pc_interactive <- function(x, ...) {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required to display an interactive plot. ",
         "Install with: install.packages('plotly')")
  p_plotly <- plotly::ggplotly(x$plot, tooltip = "text")
  print(p_plotly)
  invisible(x)
}


# ---- Internal helpers ------------------------------------------------------

## Detect which column is the criterion (HSD / LSD / Bonferroni).
#' @noRd
.pc_crit_col <- function(res) {
  cands <- intersect(names(res), c("HSD", "LSD", "Bonferroni"))
  if (length(cands) == 0L)
    stop("No criterion column (HSD, LSD, or Bonferroni) found in 'res'. ",
         "Was it produced by compare()?")
  cands[1L]
}

## Identify (a) the group columns and (b) the comparison variable.
##
## Strategy:
##   - Remove known non-factor bookkeeping columns.
##   - Count unique values per candidate column across the whole dataset.
##   - The comparison variable is the column with the most unique values
##     (typically varieties: 20–200 unique levels).
##   - All other columns with fewer unique values are treated as by-columns
##     (typically sites/treatments: 2–20 unique levels).
##
##   This cardinality approach is robust when the criterion value happens
##   to be identical across multiple groups — the criterion-based grouping
##   used previously would then fail to distinguish group vs comparison columns.
#' @noRd
.pc_detect_cols <- function(res, crit_col) {

  skip       <- c("predicted.value", "std.error", "avsed", crit_col)
  candidates <- setdiff(names(res), skip)

  if (length(candidates) == 0L)
    stop("No factor columns found in 'res' to use as the comparison variable.")

  # Count unique values for each candidate column
  n_unique <- vapply(candidates,
                     function(col) length(unique(as.character(res[[col]]))),
                     integer(1L))

  max_u    <- max(n_unique)
  comp_col <- candidates[n_unique == max_u]

  # If there is a tie for highest cardinality, warn and use the first
  if (length(comp_col) > 1L) {
    warning("Multiple columns share the highest cardinality (",
            max_u, " unique values): ",
            paste(comp_col, collapse = ", "),
            ". Using '", comp_col[1L], "' as the comparison variable.",
            call. = FALSE)
    comp_col <- comp_col[1L]
  }

  by_cols <- candidates[n_unique < max_u]

  list(by_cols  = by_cols,
       comp_col = comp_col)
}

## Build a composite group label column ("All" when no by-columns).
#' @noRd
.pc_group_label <- function(res, by_cols) {
  if (length(by_cols) == 0L)
    return(rep("All", nrow(res)))
  if (length(by_cols) == 1L)
    return(as.character(res[[by_cols]]))
  apply(res[, by_cols, drop = FALSE], 1L, paste, collapse = ":")
}

## Compact letter display (CLD) for one group.
##
## Algorithm:
##   1. Sort varieties by predicted value descending.
##   2. Start with all varieties sharing letter "a".
##   3. For each pair (i, j) with i < j (in sorted order):
##      if abs(pred_i - pred_j) > criterion, ensure they do NOT share
##      a letter: if they share all letters, assign a new letter to j
##      and all varieties between i and j that are also not significantly
##      different from j.
##   This is the classical "sweep" absorption algorithm.
#' @noRd
.pc_cld_group <- function(pred_vals, names_vals, criterion) {

  n <- length(pred_vals)
  if (n == 0L) return(character(0L))
  if (n == 1L) return(stats::setNames("a", names_vals))

  ord     <- order(pred_vals, decreasing = TRUE)
  pv_s    <- pred_vals[ord]
  nm_s    <- names_vals[ord]

  # Each variety starts with a set of letters (represented as integer set)
  letter_sets <- vector("list", n)
  for (k in seq_len(n)) letter_sets[[k]] <- 1L   # letter "a" = 1

  next_letter <- 2L

  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      if (abs(pv_s[i] - pv_s[j]) > criterion) {
        # i and j are significantly different — must not share a letter
        shared <- intersect(letter_sets[[i]], letter_sets[[j]])
        if (length(shared) > 0L) {
          # Give j (and its non-sig neighbours above it) a new letter
          new_let <- next_letter
          next_letter <- next_letter + 1L
          # Find all varieties between i+1 and j that are NOT sig diff from j
          for (k in seq(i + 1L, j)) {
            if (abs(pv_s[k] - pv_s[j]) <= criterion)
              letter_sets[[k]] <- union(letter_sets[[k]], new_let)
          }
          # Remove the shared letter from j if i still has it
          for (k in seq(i, j)) {
            if (any(letter_sets[[k]] %in% letter_sets[[i]]) &&
                abs(pv_s[k] - pv_s[j]) > criterion) {
              letter_sets[[k]] <- setdiff(letter_sets[[k]], shared)
              if (length(letter_sets[[k]]) == 0L)
                letter_sets[[k]] <- new_let
            }
          }
        }
      }
    }
  }

  # Convert integer sets to letter strings
  int_to_letters <- function(ints) {
    ints <- sort(unique(ints))
    paste(letters[ints], collapse = "")
  }
  letter_strs <- vapply(letter_sets, int_to_letters, character(1L))
  stats::setNames(letter_strs, nm_s)
}

## Build the tidy data frame for the dotplot type.
##
## When reference = NULL (default) the band is anchored to the top-ranked
## variety in each group (one-sided: band = [top - crit, top]).
## sig is a two-level logical (FALSE = not sig from best, TRUE = sig).
##
## When reference is a variety name the band is centred on that variety
## (two-sided: band = [ref - crit, ref + crit]).
## sig becomes a three-level character: "better" / "ns" / "worse".
## If reference is not found in a group a warning is issued and the
## group falls back to top-ranked behaviour.
#' @noRd
.pc_dotplot_data <- function(res, crit_col, comp_col, by_cols,
                              reference = NULL) {

  res$Group   <- .pc_group_label(res, by_cols)
  res$Variety <- as.character(res[[comp_col]])
  res$pred    <- res$predicted.value
  res$crit    <- res[[crit_col]]

  # Within each group rank by predicted value (1 = highest)
  res$rank <- ave(res$pred, res$Group,
                  FUN = function(x) rank(-x, ties.method = "first"))

  grps <- unique(res$Group)

  # Initialise output columns
  res$ref_pred <- NA_real_
  res$band_lo  <- NA_real_
  res$band_hi  <- NA_real_
  res$sig      <- NA_character_

  for (g in grps) {
    idx <- res$Group == g
    sub <- res[idx, ]
    cr  <- unique(sub$crit)

    if (!is.null(reference)) {
      ref_idx <- match(reference, sub$Variety)
      if (is.na(ref_idx)) {
        warning("Reference '", reference, "' not found in group '", g,
                "'. Falling back to top-ranked variety.", call. = FALSE)
        use_ref <- FALSE
      } else {
        use_ref <- TRUE
      }
    } else {
      use_ref <- FALSE
    }

    if (use_ref) {
      ref_pred <- sub$pred[ref_idx]
      lo       <- ref_pred - cr
      hi       <- ref_pred + cr
      sig_vec  <- ifelse(sub$pred > hi, "better",
                  ifelse(sub$pred < lo, "worse", "ns"))
    } else {
      # Default: anchor to the best
      ref_pred <- max(sub$pred, na.rm = TRUE)
      lo       <- ref_pred - cr
      hi       <- ref_pred          # one-sided
      sig_vec  <- ifelse(abs(sub$pred - ref_pred) > cr, "sig", "ns")
    }

    res$ref_pred[idx] <- ref_pred
    res$band_lo[idx]  <- lo
    res$band_hi[idx]  <- hi
    res$sig[idx]      <- sig_vec
  }

  # Backwards-compatible logical 'sig' flag used by letters/tests
  # (TRUE = any kind of significant departure from reference)
  res$sig_any <- res$sig %in% c("sig", "better", "worse")

  # Pre-compute a per-group numeric y-position (1 = lowest pred, n = highest).
  # Using a continuous numeric y-axis (with custom labels) instead of a
  # discrete factor axis means geom_rect works with finite ymin/ymax on a
  # continuous scale — critical for plotly compatibility.
  res$y_pos <- ave(res$pred, res$Group,
                   FUN = function(x) rank(x, ties.method = "first"))

  # n per group, needed for geom_rect y bounds
  res$n_group <- as.numeric(ave(seq_len(nrow(res)), res$Group,
                                FUN = length))

  res[order(res$Group, res$rank), ]
}

## Build the tidy data frame for the letters type.
#' @noRd
.pc_letters_data <- function(res, crit_col, comp_col, by_cols) {

  base <- .pc_dotplot_data(res, crit_col, comp_col, by_cols, reference = NULL)

  # Compute CLD per group
  grps <- unique(base$Group)
  letter_rows <- lapply(grps, function(g) {
    sub  <- base[base$Group == g, ]
    lmap <- .pc_cld_group(sub$pred, sub$Variety, unique(sub$crit))
    data.frame(Group   = g,
               Variety = names(lmap),
               letter  = unname(lmap),
               stringsAsFactors = FALSE)
  })
  letter_df <- do.call(rbind, letter_rows)

  merge(base, letter_df, by = c("Group", "Variety"), all.x = TRUE)
}

## Build the tidy data frame for the heatmap type.
#' @noRd
.pc_heatmap_data <- function(res, crit_col, comp_col, by_cols) {

  res$Group   <- .pc_group_label(res, by_cols)
  res$Variety <- as.character(res[[comp_col]])
  res$pred    <- res$predicted.value
  res$crit    <- res[[crit_col]]

  grps <- unique(res$Group)

  rows <- lapply(grps, function(g) {
    sub  <- res[res$Group == g, ]
    # Sort by predicted value descending for axis ordering
    sub  <- sub[order(sub$pred, decreasing = TRUE), ]
    vars <- sub$Variety
    cr   <- unique(sub$crit)
    pv   <- stats::setNames(sub$pred, vars)

    pairs <- expand.grid(Var_x = vars, Var_y = vars,
                         stringsAsFactors = FALSE)
    pairs$pred_x <- pv[pairs$Var_x]
    pairs$pred_y <- pv[pairs$Var_y]
    pairs$diff   <- abs(pairs$pred_x - pairs$pred_y)
    pairs$sig    <- pairs$diff > cr
    pairs$Group  <- g
    pairs$crit   <- cr

    # Factor levels in descending predicted value order.
    # unique() guards against duplicated levels if variety names were
    # non-unique within a group (which should not happen but is defensive).
    vars_u <- unique(vars)
    pairs$Var_x  <- factor(pairs$Var_x, levels = vars_u)
    pairs$Var_y  <- factor(pairs$Var_y, levels = rev(vars_u))
    pairs
  })

  do.call(rbind, rows)
}


# ---- Plot builders ---------------------------------------------------------

#' @noRd
.pc_plot_dotplot <- function(df, crit_col, reference, n_groups, theme, ...) {

  using_ref   <- !is.null(reference)
  has_3levels <- using_ref && any(df$sig == "better", na.rm = TRUE) ||
                              any(df$sig == "worse",  na.rm = TRUE)

  # Band data frame — uses pre-computed band_lo / band_hi.
  # ymin/ymax are finite values that span the full panel on a discrete y-axis
  # (positions run from 0.5 to n+0.5).  This avoids passing -Inf/Inf to
  # plotly, which cannot do arithmetic on infinite values.
  band_df        <- unique(df[, c("Group", "ref_pred", "band_lo",
                                  "band_hi", "n_group")])
  band_df$y_lo   <- 0.5
  band_df$y_hi   <- band_df$n_group + 0.5

  # Colour scheme and legend labels
  if (using_ref) {
    col_vals   <- c(better = "#59A14F", ns = "#4E79A7", worse = "#E15759")
    col_labels <- c(better = paste0("Sig. better than ", reference,
                                    " (", crit_col, ")"),
                    ns     = paste0("Not sig. from ", reference),
                    worse  = paste0("Sig. worse than ", reference,
                                    " (", crit_col, ")"))
    caption <- paste0("Shaded band: \u00b1", crit_col,
                      " criterion around reference variety (", reference, ").")
  } else {
    col_vals   <- c(ns = "#4E79A7", sig = "#E15759")
    col_labels <- c(ns  = paste0("Not sig. from best (", crit_col, ")"),
                    sig = paste0("Sig. from best (", crit_col, ")"))
    caption <- paste0("Shaded band: ", crit_col,
                      " criterion below the top-ranked variety. ",
                      "Red points fall outside the band.")
  }

  # Tooltip text for plotly
  df$tooltip_text <- paste0(
    "<b>", df$Variety, "</b>",
    "<br>Predicted: ", round(df$pred, 2L),
    "<br>Group: ",     df$Group,
    "<br>",            ifelse(df$sig == "ns", "Not significant",
                       ifelse(df$sig == "sig", paste0("Sig. from best (", crit_col, ")"),
                       paste0("Sig. ", df$sig, " (", crit_col, ")")))
  )

  # y-axis label lookup: numeric position -> variety name, per group
  # We build a list: group -> named vector (pos -> name)
  grp_labels <- lapply(unique(df$Group), function(g) {
    sub <- df[df$Group == g, ]
    stats::setNames(sub$Variety, sub$y_pos)
  })
  names(grp_labels) <- unique(df$Group)

  # For a single group, set up the label function for scale_y_continuous
  # (for multi-group free-scale facets each panel gets its own labels
  # automatically since y_pos 1..n is consistent within each group)
  all_y   <- sort(unique(df$y_pos))
  all_lbl <- grp_labels[[1L]][as.character(all_y)]  # use group 1 as template
  # for multi-group the varieties differ per panel — use position-based labels
  # that work across all groups (same position = same rank within group)
  label_fn <- if (n_groups == 1L) {
    function(breaks) {
      lbl <- grp_labels[[1L]][as.character(breaks)]
      ifelse(is.na(lbl), "", lbl)
    }
  } else {
    # With facet free-scales each panel has its own range; ggplot2 calls the
    # label function per panel with the panel's break values.  Since y_pos
    # maps 1..n within each group the same label function applies globally.
    function(breaks) {
      # Look up across all groups combined
      combined <- unlist(grp_labels)
      lbl <- combined[as.character(breaks)]
      # If duplicated (same pos in different groups) just take first
      ifelse(is.na(lbl), "", lbl)
    }
  }

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = pred, y = y_pos)) +
    # Criterion band on continuous y-axis — fully plotly-compatible
    ggplot2::geom_rect(
      data        = band_df,
      ggplot2::aes(xmin = band_lo, xmax = band_hi,
                   ymin = y_lo,   ymax = y_hi),
      fill        = "#FFEFC0",
      alpha       = 0.6,
      inherit.aes = FALSE
    ) +
    # Reference line
    ggplot2::geom_vline(
      data      = band_df,
      ggplot2::aes(xintercept = ref_pred),
      linetype  = "dashed",
      linewidth = 0.45,
      colour    = "grey40"
    ) +
    # All points — reference variety included so it stays at its correct
    # y_pos; the diamond overlay drawn below visually replaces the dot
    ggplot2::geom_point(
      ggplot2::aes(colour = sig, text = tooltip_text),
      size = 2.4,
      ...
    ) +
    ggplot2::scale_colour_manual(
      values = col_vals,
      labels = col_labels,
      name   = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = all_y,
      labels = label_fn,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::labs(
      x       = "Predicted value",
      y       = NULL,
      caption = caption
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom",
      panel.spacing    = ggplot2::unit(0.5, "lines"),
      axis.text.y      = ggplot2::element_text(size = 7)
    )

  # Reference variety: diamond overlay at its correct y_pos
  if (using_ref) {
    ref_rows <- df[df$Variety == reference, ]
    if (nrow(ref_rows) > 0L) {
      ref_rows$tooltip_text <- paste0(
        "<b>", ref_rows$Variety, "</b> (reference)",
        "<br>Predicted: ", round(ref_rows$pred, 2L),
        "<br>Group: ", ref_rows$Group
      )
      p <- p +
        ggplot2::geom_point(
          data        = ref_rows,
          ggplot2::aes(x = pred, y = y_pos, text = tooltip_text),
          shape       = 18L,
          size        = 4.5,
          colour      = "grey20",
          inherit.aes = FALSE
        )
    }
  }

  if (n_groups > 1L)
    p <- p + ggplot2::facet_wrap(~ Group, scales = "free")

  p
}

#' @noRd
.pc_plot_letters <- function(df, crit_col, n_groups, theme, ...) {

  df$tooltip_text <- paste0(
    "<b>", df$Variety, "</b>",
    "<br>Predicted: ", round(df$pred, 2L),
    "<br>Letter(s): ", df$letter,
    "<br>Group: ",     df$Group
  )

  # Build per-group y-axis label lookup (same approach as dotplot)
  grp_labels_l <- lapply(unique(df$Group), function(g) {
    sub <- df[df$Group == g, ]
    stats::setNames(sub$Variety, sub$y_pos)
  })
  names(grp_labels_l) <- unique(df$Group)
  all_y_l   <- sort(unique(df$y_pos))
  label_fn_l <- function(breaks) {
    combined <- unlist(grp_labels_l)
    lbl <- combined[as.character(breaks)]
    ifelse(is.na(lbl), "", lbl)
  }

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = pred, y = y_pos)) +
    ggplot2::geom_vline(
      data      = unique(df[, c("Group", "ref_pred")]),
      ggplot2::aes(xintercept = ref_pred),
      linetype  = "dashed",
      linewidth = 0.45,
      colour    = "grey60"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(text = tooltip_text),
      colour = "#4E79A7", size = 2.4, ...
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = letter),
      hjust    = -0.5,
      size     = 3.0,
      colour   = "grey20",
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      breaks = all_y_l,
      labels = label_fn_l,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::labs(
      x       = "Predicted value",
      y       = NULL,
      caption = paste0("Letters from compact letter display (",
                       crit_col, " criterion). ",
                       "Varieties sharing a letter are not significantly different.")
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing    = ggplot2::unit(0.5, "lines"),
      axis.text.y      = ggplot2::element_text(size = 7)
    )

  if (n_groups > 1L)
    p <- p + ggplot2::facet_wrap(~ Group, scales = "free")

  p
}

#' @noRd
.pc_plot_heatmap <- function(df, crit_col, n_groups, theme, ...) {

  df$tooltip_text <- paste0(
    "<b>", df$Var_x, " vs ", df$Var_y, "</b>",
    "<br>|Difference|: ", round(df$diff, 2L),
    "<br>Significant: ", ifelse(df$sig, "Yes", "No"),
    "<br>Group: ", df$Group
  )

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = Var_x, y = Var_y, fill = diff,
                      text = tooltip_text)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_point(
      data        = df[df$sig, ],
      ggplot2::aes(x = Var_x, y = Var_y),
      shape       = 4L,          # cross
      size        = 1.4,
      colour      = "white",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low  = "#E8F4FC",
      high = "#08306B",
      name = "|Difference|"
    ) +
    ggplot2::labs(
      x       = NULL,
      y       = NULL,
      caption = paste0("Tile colour: absolute pairwise difference. ",
                       "White \u00d7 marks: significant at ",
                       crit_col, " criterion.")
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1,
                                               size  = 7),
      axis.text.y      = ggplot2::element_text(size = 7),
      panel.spacing    = ggplot2::unit(0.5, "lines"),
      legend.position  = "right"
    )

  if (n_groups > 1L)
    p <- p + ggplot2::facet_wrap(~ Group, scales = "free")

  p
}


# ---- Public function -------------------------------------------------------

#' Plot Output from compare()
#'
#' @description
#' Produces one of three ggplot2 visualisations from a [compare()] result
#' object, all correctly handling the `by`-group structure (including
#' multi-factor groups) through automatic faceting.
#'
#' @details
#' The three `type` options are:
#' \describe{
#'   \item{`"dotplot"`}{Varieties sorted by predicted value (highest at top)
#'     within each group.  A shaded band of width equal to the comparison
#'     criterion is drawn below the top-ranked variety — any variety whose
#'     point falls within the band is **not** significantly different from the
#'     best.  Points outside the band are coloured red.  The dashed vertical
#'     line marks the top-ranked predicted value.}
#'   \item{`"letters"`}{Compact letter display (CLD) overlaid on the sorted
#'     dot plot.  Varieties sharing at least one letter are not significantly
#'     different.  Letters are assigned by the standard descending sweep
#'     algorithm using the comparison criterion as the significance threshold.}
#'   \item{`"heatmap"`}{An \eqn{n \times n} tile matrix with varieties sorted
#'     by predicted value on both axes.  Tile colour encodes the absolute
#'     pairwise difference; a white cross (\eqn{\times}) marks pairs that are
#'     significantly different.  Useful for large variety sets where a letter
#'     display becomes unreadable.}
#' }
#'
#' **Group detection** — `plot_compare()` infers which columns are grouping
#' variables and which is the comparison variable directly from `res`.
#' Factor columns that are constant within every group (i.e. within every
#' block of rows sharing the same criterion value) are treated as `by`
#' columns; the remaining factor column is the comparison variable.  When
#' multiple `by` columns are present their values are pasted into a
#' composite facet label.
#'
#' @param res         A data frame returned by [compare()].
#' @param type        Character string selecting the plot type. One of
#'   `"dotplot"` (default), `"letters"`, or `"heatmap"`.
#' @param reference   Character string or `NULL`. Applies to `type =
#'   "dotplot"` only. When `NULL` (default), the criterion band is anchored
#'   to the top-ranked variety in each group (one-sided band; points outside
#'   the band are red). When set to a variety name, the band is centred on
#'   that variety (`[ref \eqn{-} crit,\; ref + crit]`), the reference variety
#'   is shown as a diamond, and points are coloured green (significantly
#'   better), blue (not significant), or red (significantly worse). If the
#'   reference variety is not found in a particular group a warning is issued
#'   and that group falls back to top-ranked behaviour. Ignored with a warning
#'   for `type = "letters"` or `"heatmap"`.
#' @param interactive Logical. If `TRUE`, converts the finished ggplot to an
#'   interactive [plotly::ggplotly()] object with hover tooltips showing the
#'   variety name, predicted value, group, and significance status. Requires
#'   the \pkg{plotly} package. Default `FALSE`.
#' @param theme       A complete ggplot2 theme object.
#'   Default [ggplot2::theme_bw()].
#' @param return_data Logical. If `TRUE` returns the tidy data frame used to
#'   build the plot rather than the plot itself. Default `FALSE`.
#' @param ...         Additional arguments passed to the background
#'   `geom_point()` call (e.g. `size`, `alpha`).
#'
#' @return A `ggplot` object (when `interactive = FALSE` and
#'   `return_data = FALSE`), a `pc_interactive` object (when
#'   `interactive = TRUE` — auto-converts to plotly on `print()`; supports
#'   `+` for ggplot2 layer extensions before conversion), or a `data.frame`
#'   (when `return_data = TRUE`).
#'
#' @seealso [compare()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#' res <- compare(model,
#'                term = "Treatment:Site:Variety",
#'                by   = "Site",
#'                type = "HSD")
#'
#' # Sorted dot plot with criterion band anchored to best variety
#' plot_compare(res)
#'
#' # Anchor band to a specific check variety
#' plot_compare(res, reference = "Scout")
#'
#' # Interactive plotly version with hover tooltips
#' plot_compare(res, interactive = TRUE)
#'
#' # Compact letter display
#' plot_compare(res, type = "letters")
#'
#' # Interactive heatmap
#' plot_compare(res, type = "heatmap", interactive = TRUE)
#'
#' # Retrieve the tidy data frame
#' df <- plot_compare(res, return_data = TRUE)
#'
#' # Extend with a custom title
#' plot_compare(res, type = "letters") +
#'   ggplot2::ggtitle("Variety comparison — HSD (Tukey)")
#' }
#'
#' @export
plot_compare <- function(res,
                         type        = c("dotplot", "letters", "heatmap"),
                         reference   = NULL,
                         interactive = FALSE,
                         theme       = ggplot2::theme_bw(),
                         return_data = FALSE,
                         ...) {

  # ---- Validate -----------------------------------------------------------
  type <- match.arg(type)

  if (!is.data.frame(res))
    stop("'res' must be a data frame returned by compare().")

  required <- c("predicted.value", "std.error", "avsed")
  missing_cols <- setdiff(required, names(res))
  if (length(missing_cols))
    stop("'res' is missing expected column(s): ",
         paste(missing_cols, collapse = ", "),
         ". Was it produced by compare()?")

  if (!is.null(reference)) {
    if (!is.character(reference) || length(reference) != 1L)
      stop("'reference' must be a single character string or NULL.")
    if (type != "dotplot") {
      warning("'reference' is only used for type = \"dotplot\". Ignoring.",
              call. = FALSE)
      reference <- NULL
    }
  }

  if (!is.logical(interactive) || length(interactive) != 1L)
    stop("'interactive' must be a single logical value.")
  if (interactive && !requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required for interactive = TRUE. ",
         "Install with: install.packages('plotly')")

  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object, e.g. ggplot2::theme_bw().")

  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("'return_data' must be a single logical value.")

  # ---- Detect structure ---------------------------------------------------
  crit_col <- .pc_crit_col(res)
  cols     <- .pc_detect_cols(res, crit_col)
  by_cols  <- cols$by_cols
  comp_col <- cols$comp_col
  n_groups <- length(unique(.pc_group_label(res, by_cols)))

  # ---- Build tidy data ----------------------------------------------------
  df <- switch(type,
    dotplot = .pc_dotplot_data(res, crit_col, comp_col, by_cols, reference),
    letters = .pc_letters_data(res, crit_col, comp_col, by_cols),
    heatmap = .pc_heatmap_data(res, crit_col, comp_col, by_cols)
  )

  if (return_data) return(df)

  # ---- Build plot ---------------------------------------------------------
  # Suppress "Ignoring unknown aesthetics: text" — the text aes is only
  # consumed by plotly::ggplotly(); ggplot2 silently drops it when
  # interactive = FALSE.  We suppress the warning to keep the static
  # rendering clean.
  p <- suppressWarnings(switch(type,
    dotplot = .pc_plot_dotplot(df, crit_col, reference, n_groups, theme, ...),
    letters = .pc_plot_letters(df, crit_col, n_groups, theme, ...),
    heatmap = .pc_plot_heatmap(df, crit_col, n_groups, theme, ...)
  ))

  # ---- Optionally wrap for lazy plotly conversion -------------------------
  # Return a pc_interactive wrapper so the user can still use + to add
  # ggplot2 layers before the plotly conversion happens at print() time.
  if (interactive)
    return(.pc_make_interactive(p))

  p
}
