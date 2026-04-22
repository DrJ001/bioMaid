# ============================================================
#  plot_randomRegress.R
#  ggplot2 visualisation for randomRegress() output.
# ============================================================

# ---- Data preparation helpers ------------------------------------------

#' @noRd
.rreg_regress_data <- function(res, treatments) {

  blups       <- res$blups
  cond_list   <- res$cond_list
  conditioned <- names(Filter(Negate(is.null), cond_list))

  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  rows <- lapply(conditioned, function(lv_j) {
    A_j     <- cond_list[[lv_j]]
    cond_lv <- A_j[1L]
    beta_j  <- res$beta[[lv_j]]          # ns x |A_j|, rownames = sites
    site_beta <- setNames(beta_j[, cond_lv], rownames(beta_j))
    data.frame(
      Site       = blups$Site,
      Variety    = blups$Variety,
      x          = blups[[cond_lv]],
      y          = blups[[lv_j]],
      pair_label = paste0(lv_j, " | ", cond_lv),
      beta       = site_beta[as.character(blups$Site)],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @noRd
.rreg_quadrant_data <- function(res, treatments) {

  blups       <- res$blups
  cond_list   <- res$cond_list
  conditioned <- names(Filter(Negate(is.null), cond_list))

  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  eff_lv <- names(Filter(is.null, cond_list))[1L]

  rows <- lapply(conditioned, function(lv_j) {
    resp_col <- paste0("resp.", lv_j)
    if (!(resp_col %in% names(blups)))
      stop("Column '", resp_col, "' not found in res$blups.")
    cond_lv <- cond_list[[lv_j]][1L]
    data.frame(
      Site       = blups$Site,
      Variety    = blups$Variety,
      x          = blups[[eff_lv]],
      y          = blups[[resp_col]],
      pair_label = paste0(lv_j, " | ", cond_lv),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @noRd
.rreg_gmat_data <- function(res) {

  Gmat <- res$Gmat
  if (is.null(Gmat)) stop("res$Gmat is NULL.")

  cormat <- cov2cor(Gmat)
  nms    <- rownames(cormat)
  n      <- length(nms)

  rows <- vector("list", n * n)
  idx  <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      rows[[idx]] <- data.frame(row_var = nms[i], col_var = nms[j],
                                corr = cormat[i, j],
                                stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  df         <- do.call(rbind, rows)
  df$row_var <- factor(df$row_var, levels = nms)
  df$col_var <- factor(df$col_var, levels = rev(nms))
  df
}

# ---- Default highlight selection ---------------------------------------

# Selects 6 varieties to annotate: 3 from the top-right quadrant and 3 from
# the bottom-left quadrant of the efficiency x responsiveness space.
# Within each quadrant the selection uses polar angles to represent three
# archetypes: efficiency (angle ~ 0), balanced (angle ~ pi/4), and
# responsiveness (angle ~ pi/2).
# @param qdata  Data frame from .rreg_quadrant_data()
# @return Data frame with columns Variety and group ("tr" or "bl")
#' @noRd
.rreg_default_highlights <- function(qdata) {

  # Collapse to one representative point per variety
  var_means <- aggregate(cbind(x, y) ~ Variety, data = qdata, FUN = mean)

  .pick_quad <- function(df, sx, sy, grp) {
    in_q <- df[df$x * sx > 0 & df$y * sy > 0, , drop = FALSE]
    if (nrow(in_q) == 0L)
      return(data.frame(Variety = character(0L), group = character(0L),
                        stringsAsFactors = FALSE))

    in_q$dist  <- sqrt(in_q$x^2 + in_q$y^2)
    in_q$angle <- atan2(abs(in_q$y), abs(in_q$x))

    # Filter to > 50th percentile distance; relax progressively if needed
    for (pct in c(0.5, 0.33, 0.0)) {
      cands <- in_q[in_q$dist > quantile(in_q$dist, pct), , drop = FALSE]
      if (nrow(cands) >= 3L) break
    }

    # Pick three archetypes sequentially, removing each pick from the pool
    selectors <- list(
      function(d) d$Variety[which.min(d$angle)],              # efficiency
      function(d) d$Variety[which.min(abs(d$angle - pi/4))],  # balanced
      function(d) d$Variety[which.max(d$angle)]               # responsiveness
    )

    chosen <- character(0L)
    pool   <- cands
    for (sel in selectors) {
      if (nrow(pool) == 0L) break
      pick   <- sel(pool)
      chosen <- c(chosen, pick)
      pool   <- pool[pool$Variety != pick, , drop = FALSE]
    }

    data.frame(Variety = chosen, group = grp, stringsAsFactors = FALSE)
  }

  rbind(
    .pick_quad(var_means,  1,  1, "tr"),
    .pick_quad(var_means, -1, -1, "bl")
  )
}

# ---- Shared highlight layer builder ------------------------------------

# Returns a list of ggplot2 layers (geom_point + geom_text) for highlighted
# varieties.  Returns an empty list when hl is NULL.
#' @noRd
.rreg_highlight_layers <- function(df, hl, label_col = "Variety") {

  if (is.null(hl) || nrow(hl) == 0L) return(list())

  hl_data <- merge(df, hl, by = "Variety")

  # Colour palette: top-right warm, bottom-left cool, custom green
  col_map <- c(tr = "#E07B39", bl = "#4E79A7", custom = "#59A14F")

  list(
    ggplot2::geom_point(
      data        = hl_data[hl_data$group == "tr", ],
      ggplot2::aes(x = x, y = y),
      colour      = col_map["tr"],
      size        = 2.8,
      inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data        = hl_data[hl_data$group == "bl", ],
      ggplot2::aes(x = x, y = y),
      colour      = col_map["bl"],
      size        = 2.8,
      inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data        = hl_data[hl_data$group == "custom", ],
      ggplot2::aes(x = x, y = y),
      colour      = col_map["custom"],
      size        = 2.8,
      inherit.aes = FALSE
    ),
    ggplot2::geom_text(
      data         = hl_data,
      ggplot2::aes(x = x, y = y, label = Variety,
                   colour = group),
      size         = 2.5,
      vjust        = -0.7,
      fontface     = "bold",
      inherit.aes  = FALSE,
      show.legend  = FALSE
    ),
    ggplot2::scale_colour_manual(
      values = col_map,
      guide  = "none"
    )
  )
}

# ---- Plot builders -------------------------------------------------------

#' @noRd
.rreg_plot_regress <- function(df, hl, theme, ...) {

  panel_df <- unique(df[, c("pair_label", "Site", "beta")])
  panel_df$beta_label <- paste0("\u03b2 = ", round(panel_df$beta, 2L))

  base_col <- if (is.null(hl)) "#4E79A7" else "grey78"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_point(size = 1.6, alpha = 0.7, colour = base_col, ...) +
    .rreg_highlight_layers(df, hl) +
    ggplot2::geom_abline(
      data     = panel_df,
      ggplot2::aes(slope = beta, intercept = 0),
      linetype = "dotted", linewidth = 0.7, colour = "grey20"
    ) +
    ggplot2::geom_text(
      data        = panel_df,
      ggplot2::aes(x = -Inf, y = Inf, label = beta_label),
      hjust       = -0.15, vjust = 1.5,
      size        = 3.2, colour = "grey20", fontface = "italic",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(pair_label ~ Site, scales = "free") +
    ggplot2::labs(
      x       = "Conditioning treatment BLUP",
      y       = "Conditioned treatment BLUP",
      caption = "Dotted line: site-specific random regression"
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing    = ggplot2::unit(0.4, "lines")
    )
  p
}

#' @noRd
.rreg_plot_quadrant <- function(df, hl, theme, ...) {

  base_col <- if (is.null(hl)) "#E15759" else "grey78"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    ggplot2::geom_point(size = 1.6, alpha = 0.7, colour = base_col, ...) +
    .rreg_highlight_layers(df, hl) +
    ggplot2::facet_grid(pair_label ~ Site, scales = "free") +
    ggplot2::labs(
      x       = "Efficiency (unconditional BLUP)",
      y       = "Responsiveness (conditional BLUP)",
      caption = "Dotted lines at zero divide each panel into four quadrants"
    ) +
    theme +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.spacing    = ggplot2::unit(0.4, "lines")
    )
  p
}

#' @noRd
.rreg_plot_gmat <- function(df, theme, ...) {

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_var, y = row_var,
                                         fill = corr)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(abs(corr) > 0.01, round(corr, 2L), "")),
      size = 2.5, colour = "grey20"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575B4", mid = "white", high = "#D73027",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7)
    )
  p
}


# ---- Public function -----------------------------------------------------

#' Plot Output from randomRegress()
#'
#' @description
#' Produces one of three ggplot2 visualisations from a [randomRegress()] result
#' object.  All plot types return a `ggplot` object that can be further
#' customised with the standard `+` operator.
#'
#' @details
#' The three `type` options are:
#' \describe{
#'   \item{`"regress"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     site (columns).  Each panel plots the raw conditioned-treatment BLUPs
#'     (y) against the conditioning-treatment BLUPs (x) for one site x one
#'     treatment pair.  A dotted random regression line with the site-specific
#'     beta slope passes through the origin.  The site-specific
#'     \eqn{\hat{\beta}} is annotated in the top-left corner of each panel.}
#'   \item{`"quadrant"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     site (columns).  Each panel plots responsiveness BLUPs (y) against the
#'     conditioning treatment BLUPs (x = efficiency) for one site x one
#'     treatment pair.  Dotted zero reference lines on both axes divide each
#'     panel into four quadrants.}
#'   \item{`"gmat"`}{Heatmap of the G-matrix converted to a correlation
#'     matrix via [stats::cov2cor()].  Fill uses a diverging palette centred
#'     at zero.}
#' }
#'
#' **Variety highlighting** (`type = "regress"` and `"quadrant"` only):
#' By default, six varieties are identified and annotated across all panels.
#' Three come from the top-right quadrant of the efficiency x responsiveness
#' space (above average on both axes) and three from the bottom-left quadrant
#' (below average on both axes).  Within each quadrant the three varieties are
#' chosen by their polar angle to represent an efficiency archetype (small
#' angle), a balanced archetype (angle near \eqn{\pi/4}), and a
#' responsiveness archetype (large angle), selecting only from varieties whose
#' distance from the origin exceeds the 50th percentile within their quadrant.
#' The same six varieties are consistently annotated across all site and
#' treatment-pair panels.
#'
#' @param res         A list returned by [randomRegress()].
#' @param type        Character string selecting the plot type. One of
#'   `"regress"` (default), `"quadrant"`, or `"gmat"`.
#' @param treatments  Character vector restricting which conditioned treatments
#'   are included in the plot.  `NULL` (default) includes all conditioned
#'   treatments.  Ignored for `type = "gmat"`.
#' @param highlight   Controls variety annotation for `"regress"` and
#'   `"quadrant"` plots. One of:
#'   \describe{
#'     \item{`"default"`}{Automatically selects 6 varieties (3 top-right,
#'       3 bottom-left) using the polar-angle algorithm described in Details.}
#'     \item{character vector}{Highlight exactly the named varieties.}
#'     \item{`NULL`}{No highlighting; all points drawn in a single colour.}
#'   }
#'   Default is `"default"`.  Ignored for `type = "gmat"`.
#' @param theme       A complete ggplot2 theme object. Default
#'   [ggplot2::theme_bw()].
#' @param return_data Logical. If `TRUE` returns the tidy data frame used to
#'   build the plot rather than the plot itself. Default `FALSE`.
#' @param ...         Additional arguments passed to the background
#'   `geom_point()` call (e.g. `size`, `alpha`, `shape`). Not used for
#'   `"gmat"`.
#'
#' @return A `ggplot` object (when `return_data = FALSE`) or a `data.frame`
#'   (when `return_data = TRUE`).
#'
#' @seealso [randomRegress()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#' res <- randomRegress(model, levs = c("N0", "N1", "N2"))
#'
#' # Regression plot with default 6 highlighted varieties
#' plot_randomRegress(res)
#'
#' # Quadrant plot with user-specified highlights
#' plot_randomRegress(res, type = "quadrant",
#'                   highlight = c("Var01", "Var15", "Var22"))
#'
#' # Suppress highlighting
#' plot_randomRegress(res, type = "quadrant", highlight = NULL)
#'
#' # G-matrix correlation heatmap
#' plot_randomRegress(res, type = "gmat")
#'
#' # Retrieve the tidy data frame
#' df <- plot_randomRegress(res, type = "quadrant", return_data = TRUE)
#' }
#'
#' @export
plot_randomRegress <- function(res,
                               type        = c("regress", "quadrant", "gmat"),
                               treatments  = NULL,
                               highlight   = "default",
                               theme       = ggplot2::theme_bw(),
                               return_data = FALSE,
                               ...) {

  # ---- Validate ----------------------------------------------------------
  type <- match.arg(type)

  required <- c("blups", "beta", "Gmat", "cond_list", "type")
  missing  <- setdiff(required, names(res))
  if (length(missing))
    stop("'res' is missing expected element(s): ",
         paste(missing, collapse = ", "),
         ". Was it produced by randomRegress()?")

  if (!is.null(treatments) && !is.character(treatments))
    stop("'treatments' must be a character vector or NULL.")

  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object, e.g. ggplot2::theme_bw().")

  # ---- Build tidy data ---------------------------------------------------
  df <- switch(type,
    regress  = .rreg_regress_data( res, treatments),
    quadrant = .rreg_quadrant_data(res, treatments),
    gmat     = .rreg_gmat_data(    res)
  )

  if (return_data) return(df)

  # ---- Resolve highlights (regress + quadrant only) ----------------------
  hl <- NULL
  if (type %in% c("regress", "quadrant")) {
    if (identical(highlight, "default")) {
      # Always derive highlights from quadrant space
      qdata <- if (type == "quadrant") df else
                 .rreg_quadrant_data(res, treatments)
      hl    <- .rreg_default_highlights(qdata)
    } else if (is.character(highlight) && length(highlight) > 0L) {
      all_vars <- unique(as.character(df$Variety))
      bad      <- setdiff(highlight, all_vars)
      if (length(bad))
        warning("'highlight' varieties not found in data: ",
                paste(bad, collapse = ", "))
      valid <- intersect(highlight, all_vars)
      if (length(valid) > 0L)
        hl <- data.frame(Variety = valid, group = "custom",
                         stringsAsFactors = FALSE)
    }
    # NULL highlight: hl stays NULL
  }

  # ---- Build plot --------------------------------------------------------
  switch(type,
    regress  = .rreg_plot_regress( df, hl, theme, ...),
    quadrant = .rreg_plot_quadrant(df, hl, theme, ...),
    gmat     = .rreg_plot_gmat(    df,     theme, ...)
  )
}
