# ============================================================
#  plot_padTrial.R
#  Before/after field-layout tile map for padTrial() output.
#
#  Design:
#   - Two panels stacked vertically: Before (top) / After (bottom).
#   - Tiles coloured by plot type (type_col).
#   - Padded cells (add == "new") shown in light grey so they
#     stand out clearly against the original plot types.
#   - When split is supplied, panels become rows of a facet grid
#     with groups as columns: Before on top, After on bottom.
#   - data = NULL (default): Before panel reconstructed from
#     add == "old" rows — gaps appear where missing cells were.
#   - data supplied: Before panel shows the full original layout
#     including guard rows and cells outside the bounding box.
#   - Optional tile labels from any column in result / data.
#
#  Pattern mirrors other plot_* functions:
#    public function -> .pt_build_data() -> .pt_build_plot()
# ============================================================


# ---- Private helpers -------------------------------------------------------

## Parse "Row:Column" pattern string into a named list.
#' @noRd
.pt_parse_pattern <- function(pattern) {
  parts <- strsplit(pattern, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 2L)
    stop("'pattern' must be two column names separated by ':', ",
         "e.g. \"Row:Column\".")
  list(row_col = trimws(parts[1L]), col_col = trimws(parts[2L]))
}

## Build the tidy data frame combining Before and After panels.
#' @noRd
.pt_build_data <- function(result, data, type_col, row_col, col_col,
                            split, label) {

  # ---- After panel (always from result) -----------------------------------
  aft            <- result
  aft$fill_group <- ifelse(aft$add == "new",
                           "Padded",
                           as.character(aft[[type_col]]))
  aft$row_num    <- as.numeric(as.character(aft[[row_col]]))
  aft$col_num    <- as.numeric(as.character(aft[[col_col]]))
  aft$panel      <- "After"
  aft$is_new     <- aft$add == "new"

  # ---- Before panel -------------------------------------------------------
  if (!is.null(data)) {
    # Full original data supplied: show complete layout including guards
    bef            <- data
    bef$fill_group <- as.character(bef[[type_col]])
    bef$is_new     <- FALSE
  } else {
    # Reconstruct from result: only "old" rows are kept.
    # "new" rows are intentionally excluded so their positions appear as
    # white gaps in the tile map — visually showing where the holes were.
    bef            <- result[result$add == "old", , drop = FALSE]
    bef$fill_group <- as.character(bef[[type_col]])
    bef$is_new     <- FALSE
  }
  bef$row_num <- as.numeric(as.character(bef[[row_col]]))
  bef$col_num <- as.numeric(as.character(bef[[col_col]]))
  bef$panel   <- "Before"

  # ---- Group column -------------------------------------------------------
  mk_grp <- function(df) {
    if (is.null(split)) return(rep("All", nrow(df)))
    if (length(split) == 1L)
      as.character(df[[split]])
    else
      apply(df[, split, drop = FALSE], 1L, paste, collapse = " \u00b7 ")
  }
  aft$group <- mk_grp(aft)
  bef$group <- mk_grp(bef)

  # ---- Label column -------------------------------------------------------
  mk_label <- function(df, is_aft) {
    if (is.null(label)) return(rep("", nrow(df)))
    if (!label %in% names(df)) return(rep("", nrow(df)))
    txt <- as.character(df[[label]])
    # Suppress label on padded tiles in After panel
    if (is_aft) txt[df$add == "new"] <- ""
    txt
  }
  aft$label_text <- mk_label(aft, TRUE)
  bef$label_text <- mk_label(bef, FALSE)

  # ---- Combine and factorise ----------------------------------------------
  keep <- c("row_num", "col_num", "fill_group", "panel",
            "group", "is_new", "label_text")
  out       <- rbind(aft[, keep, drop = FALSE],
                     bef[, keep, drop = FALSE])
  out$panel <- factor(out$panel, levels = c("Before", "After"))
  out
}

## Build the ggplot tile map object.
#' @noRd
.pt_build_plot <- function(data, type_col, split, label, data_null,
                            theme, ...) {

  has_grps <- length(unique(data$group)) > 1L

  # ---- Colour palette -----------------------------------------------------
  # Distinct qualitative colours for observed plot types; Padded = light grey.
  type_vals <- sort(setdiff(unique(data$fill_group), "Padded"))
  qual_pal  <- c("#4E79A7", "#59A14F", "#F28E2B", "#E15759",
                 "#B07AA1", "#76B7B2", "#FF9DA7", "#9C755F",
                 "#EDC948", "#BAB0AC")
  n_types   <- length(type_vals)
  type_cols <- stats::setNames(qual_pal[seq_len(n_types)], type_vals)
  fill_pal  <- c(type_cols, Padded = "#DEDEDE")

  # ---- Caption ------------------------------------------------------------
  cap <- if (data_null)
    paste0("Before panel reconstructed from result ",
           "(original data not supplied \u2014 guard rows and excluded ",
           "cells not shown).")
  else
    NULL

  # ---- Base plot ----------------------------------------------------------
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = col_num, y = row_num, fill = fill_group)
  ) +
    # Main tile layer
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +

    # ---- Scales -----------------------------------------------------------
    ggplot2::scale_fill_manual(
      values       = fill_pal,
      name         = type_col,
      na.translate = FALSE
    ) +

    # Row 1 at the top (north), matching field convention
    ggplot2::scale_y_reverse(
      breaks = function(x) unique(as.integer(pretty(x)))
    ) +

    ggplot2::scale_x_continuous(
      breaks = function(x) unique(as.integer(pretty(x)))
    ) +

    ggplot2::labs(
      x       = "Column",
      y       = "Row",
      caption = cap
    ) +

    theme +

    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 8),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      legend.position  = "right",
      plot.caption     = ggplot2::element_text(size    = 7,
                                               colour  = "grey50",
                                               hjust   = 0)
    )

  # ---- Optional tile labels -----------------------------------------------
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label_text),
      size        = 2.2,
      colour      = "grey15",
      show.legend = FALSE,
      na.rm       = TRUE
    )
  }

  # ---- Faceting -----------------------------------------------------------
  if (has_grps) {
    # Groups as columns, Before/After as rows → Before always on top
    p <- p + ggplot2::facet_grid(
      panel ~ group,
      scales = "free",
      space  = "free"
    )
  } else {
    # No grouping: two panels stacked, Before on top
    p <- p + ggplot2::facet_wrap(
      ~ panel,
      ncol   = 1L,
      scales = "free"
    )
  }

  p
}


# ---- Main exported function ------------------------------------------------

#' Before/After Field Layout Plot for padTrial() Results
#'
#' @description
#' Produces a tile map showing the field trial layout **before** and **after**
#' [padTrial()] is applied.
#'
#' The **Before** panel shows the original layout. If `data` is supplied this
#' includes guard rows and any cells outside the bounding box; if `data` is
#' `NULL` (default) the Before panel is reconstructed from the retained
#' (`"old"`) rows of `result` and gaps appear as white space where missing
#' cells were located.
#'
#' The **After** panel shows the padded result. Inserted cells (`add == "new"`)
#' are shown in light grey so they are immediately distinguishable from
#' original plots.
#'
#' Tiles are coloured by the plot-type column. When `split` is supplied,
#' groups appear as columns in a facet grid with Before/After as rows, keeping
#' Before always on top.
#'
#' @param result      Data frame returned by [padTrial()]. Must contain an
#'   `add` column with values `"old"` and `"new"`.
#' @param data        Optional original data frame passed to [padTrial()].
#'   When supplied the Before panel shows the complete original layout
#'   including guard rows. When `NULL` (default) the Before panel is
#'   reconstructed from the `"old"` rows of `result`.
#' @param type_col    Character. Name of the plot-type column used to colour
#'   tiles. Default `"Type"`.
#' @param pattern     Character. Colon-separated names of the row and column
#'   coordinate columns, e.g. `"Row:Column"`. Default `"Row:Column"`.
#' @param split       Character vector of column name(s) matching the `split`
#'   argument used in [padTrial()]. Groups appear as columns in the facet
#'   grid. Default `NULL` (no grouping).
#' @param label       Character or `NULL`. Name of a column whose values are
#'   printed as text inside each tile. `NULL` (default) produces no labels.
#' @param theme       A ggplot2 theme object. Default [ggplot2::theme_bw()].
#' @param return_data Logical. If `TRUE` returns the tidy data frame used to
#'   build the plot. Default `FALSE`.
#' @param ...         Reserved for future use.
#'
#' @return A [ggplot2::ggplot] object returned invisibly (extendable with `+`),
#'   or a data frame when `return_data = TRUE`.
#'
#' @examples
#' \dontrun{
#' ## Basic use — Before panel reconstructed from result
#' res <- padTrial(trial_data, split = "Block", verbose = FALSE)
#' plot_padTrial(res)
#'
#' ## Supply original data to include guard rows in Before panel
#' plot_padTrial(res, data = trial_data)
#'
#' ## Multi-block split — groups as columns
#' plot_padTrial(res, data = trial_data, split = "Block")
#'
#' ## Label each tile with the genotype name
#' plot_padTrial(res, label = "Geno")
#'
#' ## Return tidy data for bespoke customisation
#' df <- plot_padTrial(res, return_data = TRUE)
#' }
#'
#' @seealso [padTrial()]
#' @export
plot_padTrial <- function(result,
                           data        = NULL,
                           type_col    = "Type",
                           pattern     = "Row:Column",
                           split       = NULL,
                           label       = NULL,
                           theme       = ggplot2::theme_bw(),
                           return_data = FALSE,
                           ...) {

  # ---- Validate result ----------------------------------------------------
  if (!is.data.frame(result))
    stop("'result' must be a data frame returned by padTrial().")
  if (!"add" %in% names(result))
    stop("'result' must contain an 'add' column. ",
         "Was it returned by padTrial()?")
  if (!all(result$add %in% c("old", "new")))
    stop("'result$add' must contain only \"old\" and \"new\".")

  # ---- Parse and validate pattern -----------------------------------------
  rc      <- .pt_parse_pattern(pattern)
  row_col <- rc$row_col
  col_col <- rc$col_col

  bad_rc <- setdiff(c(row_col, col_col), names(result))
  if (length(bad_rc))
    stop("Column(s) from 'pattern' not found in 'result': ",
         paste(bad_rc, collapse = ", "), ".")

  # ---- Validate type_col --------------------------------------------------
  if (!is.character(type_col) || length(type_col) != 1L)
    stop("'type_col' must be a single character string.")
  if (!type_col %in% names(result))
    stop("'type_col' column \"", type_col, "\" not found in 'result'.")

  # ---- Validate split -----------------------------------------------------
  if (!is.null(split)) {
    bad_sp <- setdiff(split, names(result))
    if (length(bad_sp))
      stop("'split' column(s) not found in 'result': ",
           paste(bad_sp, collapse = ", "), ".")
    if (!is.null(data)) {
      bad_sp_d <- setdiff(split, names(data))
      if (length(bad_sp_d))
        stop("'split' column(s) not found in 'data': ",
             paste(bad_sp_d, collapse = ", "), ".")
    }
  }

  # ---- Validate label -----------------------------------------------------
  if (!is.null(label)) {
    if (!is.character(label) || length(label) != 1L)
      stop("'label' must be a single character string or NULL.")
    if (!label %in% names(result))
      stop("'label' column \"", label, "\" not found in 'result'.")
  }

  # ---- Validate data ------------------------------------------------------
  if (!is.null(data)) {
    if (!is.data.frame(data))
      stop("'data' must be a data frame or NULL.")
    bad_d <- setdiff(c(row_col, col_col, type_col), names(data))
    if (length(bad_d))
      stop("Column(s) not found in 'data': ",
           paste(bad_d, collapse = ", "), ".")
  }

  # ---- Validate theme / return_data ---------------------------------------
  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object.")
  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("'return_data' must be TRUE or FALSE.")

  # ---- Build tidy data frame ----------------------------------------------
  tidy <- .pt_build_data(result,
                         data     = data,
                         type_col = type_col,
                         row_col  = row_col,
                         col_col  = col_col,
                         split    = split,
                         label    = label)

  if (return_data) return(tidy)

  # ---- Build and return plot ----------------------------------------------
  p <- .pt_build_plot(tidy,
                      type_col  = type_col,
                      split     = split,
                      label     = label,
                      data_null = is.null(data),
                      theme     = theme,
                      ...)
  invisible(p)
}
