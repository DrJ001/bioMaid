# ============================================================
#  plot_fixedRegress.R
#  ggplot2 visualisation for fixedRegress() output.
#
#  Analogous to plot_randomRegress() but for BLUEs / OLS.
#  Key differences from the random-effects version:
#   - Group column name is dynamic (inferred from blues)
#   - Genotype column name is dynamic (inferred from blues)
#   - Regression line has a non-zero OLS intercept
#   - No G-matrix; only two plot types: "regress" and "quadrant"
#  Output data frames standardise group -> "Group", genotype -> "Genotype"
#  to simplify downstream ggplot2 calls.
# ============================================================

# ---- Column detection --------------------------------------------------

# Infer the group and genotype column names from the blues data frame.
# blues col 1 = group (by_col), col 2 = genotype (gnam).
#' @noRd
.freg_detect_cols <- function(res) {
  nms <- names(res$blues)
  list(grp_col  = nms[1L],
       geno_col = nms[2L])
}

# ---- Data preparation helpers ------------------------------------------

#' @noRd
.freg_regress_data <- function(res, treatments, centre) {

  blues     <- res$blues
  cond_list <- res$cond_list
  cols      <- .freg_detect_cols(res)
  grp_col   <- cols$grp_col
  geno_col  <- cols$geno_col

  conditioned <- names(Filter(Negate(is.null), cond_list))
  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  rows <- lapply(conditioned, function(lv_j) {
    A_j     <- cond_list[[lv_j]]
    cond_lv <- A_j[1L]

    # Extract group-specific beta from the beta data frame
    beta_df   <- res$beta[[lv_j]]              # data frame: grp_col + A_j cols
    site_beta <- setNames(beta_df[[cond_lv]],
                          as.character(beta_df[[grp_col]]))

    df <- data.frame(
      Group      = as.character(blues[[grp_col]]),
      Genotype   = as.character(blues[[geno_col]]),
      x          = blues[[cond_lv]],
      y          = blues[[lv_j]],
      pair_label = paste0(lv_j, " | ", cond_lv),
      beta       = site_beta[as.character(blues[[grp_col]])],
      stringsAsFactors = FALSE
    )

    # Compute OLS intercept per group: mean(y) - beta * mean(x)
    grp_stats <- do.call(rbind, lapply(unique(df$Group), function(g) {
      sub <- df[df$Group == g, ]
      b   <- unique(sub$beta)
      data.frame(
        Group      = g,
        pair_label = unique(sub$pair_label),
        intercept  = mean(sub$y, na.rm = TRUE) - b * mean(sub$x, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))

    out <- merge(df, grp_stats[, c("Group", "pair_label", "intercept")],
                 by = c("Group", "pair_label"), all.x = TRUE)

    if (centre) {
      # Subtract within-group means from both axes.
      # beta (slope) is invariant to location shifts; intercept becomes 0.
      out$x         <- out$x - ave(out$x, out$Group, FUN = mean)
      out$y         <- out$y - ave(out$y, out$Group, FUN = mean)
      out$intercept <- 0
    }
    out
  })

  do.call(rbind, rows)
}

#' @noRd
.freg_quadrant_data <- function(res, treatments, centre) {

  blues     <- res$blues
  cond_list <- res$cond_list
  cols      <- .freg_detect_cols(res)
  grp_col   <- cols$grp_col
  geno_col  <- cols$geno_col

  conditioned <- names(Filter(Negate(is.null), cond_list))
  if (!is.null(treatments))
    conditioned <- intersect(conditioned, treatments)
  if (length(conditioned) == 0L)
    stop("No conditioned treatments to plot. Check 'treatments'.")

  eff_lv <- names(Filter(is.null, cond_list))[1L]

  rows <- lapply(conditioned, function(lv_j) {
    resp_col <- paste0("resp.", lv_j)
    if (!(resp_col %in% names(blues)))
      stop("Column '", resp_col, "' not found in res$blues.")
    cond_lv <- cond_list[[lv_j]][1L]

    df <- data.frame(
      Group      = as.character(blues[[grp_col]]),
      Genotype   = as.character(blues[[geno_col]]),
      x          = blues[[eff_lv]],
      y          = blues[[resp_col]],
      pair_label = paste0(lv_j, " | ", cond_lv),
      stringsAsFactors = FALSE
    )

    if (centre)
      df$x <- df$x - ave(df$x, df$Group, FUN = mean)

    df
  })

  do.call(rbind, rows)
}

# ---- Default highlight selection ---------------------------------------

# Same polar-angle algorithm as .rreg_default_highlights() but operating on
# the fixed-effects quadrant data frame (which uses "Genotype" not "Variety").
#' @noRd
.freg_default_highlights <- function(qdata) {

  var_means <- aggregate(cbind(x, y) ~ Genotype, data = qdata, FUN = mean)

  .pick_quad <- function(df, sx, sy, grp) {
    in_q <- df[df$x * sx > 0 & df$y * sy > 0, , drop = FALSE]
    if (nrow(in_q) == 0L)
      return(data.frame(Genotype = character(0L), group = character(0L),
                        stringsAsFactors = FALSE))

    in_q$dist  <- sqrt(in_q$x^2 + in_q$y^2)
    in_q$angle <- atan2(abs(in_q$y), abs(in_q$x))

    for (pct in c(0.5, 0.33, 0.0)) {
      cands <- in_q[in_q$dist > quantile(in_q$dist, pct), , drop = FALSE]
      if (nrow(cands) >= 3L) break
    }

    selectors <- list(
      function(d) d$Genotype[which.min(d$angle)],
      function(d) d$Genotype[which.min(abs(d$angle - pi / 4))],
      function(d) d$Genotype[which.max(d$angle)]
    )

    chosen <- character(0L)
    pool   <- cands
    for (sel in selectors) {
      if (nrow(pool) == 0L) break
      pick   <- sel(pool)
      chosen <- c(chosen, pick)
      pool   <- pool[pool$Genotype != pick, , drop = FALSE]
    }

    if (length(chosen) == 0L)
      return(data.frame(Genotype = character(0L), group = character(0L),
                        stringsAsFactors = FALSE))

    data.frame(Genotype = chosen, group = grp, stringsAsFactors = FALSE)
  }

  rbind(
    .pick_quad(var_means,  1,  1, "tr"),
    .pick_quad(var_means, -1, -1, "bl")
  )
}

# ---- Shared highlight layer builder ------------------------------------

# Analogous to .rreg_highlight_layers() but uses "Genotype" for labels.
#' @noRd
.freg_highlight_layers <- function(df, hl) {

  if (is.null(hl) || nrow(hl) == 0L) return(list())

  hl_data <- merge(df, hl, by = "Genotype")

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
      ggplot2::aes(x = x, y = y, label = Genotype, colour = group),
      size         = 2.5,
      vjust        = -0.7,
      fontface     = "bold",
      inherit.aes  = FALSE,
      show.legend  = FALSE
    ),
    ggplot2::scale_colour_manual(values = col_map, guide = "none")
  )
}

# ---- Plot builders -------------------------------------------------------

#' @noRd
.freg_plot_regress <- function(df, hl, centre, theme, ...) {

  # One row per Group x pair_label panel for ablines and beta annotation
  panel_df <- unique(df[, c("pair_label", "Group", "beta", "intercept")])
  panel_df$beta_label <- paste0("\u03b2 = ", round(panel_df$beta, 2L))

  base_col <- if (is.null(hl)) "#4E79A7" else "grey78"

  # Reference lines: zero when centred; treatment means when not centred.
  # When centre = FALSE, x and y are on absolute BLUE scale (~4500) so 0 is
  # off-screen.  Instead, draw lines at the within-group means of x and y --
  # these fall in the middle of each panel and give the same above/below
  # average split as the zero line does in the centred case.
  ref_layers <- if (centre) {
    list(
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70"),
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70")
    )
  } else {
    x_means <- aggregate(x ~ Group,              data = df, FUN = mean)
    y_means <- aggregate(y ~ Group + pair_label, data = df, FUN = mean)
    list(
      ggplot2::geom_vline(
        data      = x_means,
        ggplot2::aes(xintercept = x),
        linewidth = 0.25, colour = "grey70"
      ),
      ggplot2::geom_hline(
        data      = y_means,
        ggplot2::aes(yintercept = y),
        linewidth = 0.25, colour = "grey70"
      )
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ref_layers +
    ggplot2::geom_point(size = 1.6, alpha = 0.7, colour = base_col, ...) +
    .freg_highlight_layers(df, hl) +
    ggplot2::geom_abline(
      data     = panel_df,
      ggplot2::aes(slope = beta, intercept = intercept),
      linetype = "dotted", linewidth = 0.7, colour = "grey20"
    ) +
    ggplot2::geom_text(
      data        = panel_df,
      ggplot2::aes(x = -Inf, y = Inf, label = beta_label),
      hjust = -0.15, vjust = 1.5,
      size = 3.2, colour = "grey20", fontface = "italic",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(pair_label ~ Group, scales = "free") +
    ggplot2::labs(
      x       = if (centre) "Conditioning BLUE (centred)" else
                             "Conditioning treatment BLUE",
      y       = if (centre) "Conditioned BLUE (centred)"  else
                             "Conditioned treatment BLUE",
      caption = "Dotted line: OLS regression within group"
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
.freg_plot_quadrant <- function(df, hl, centre, theme, ...) {

  base_col <- if (is.null(hl)) "#E15759" else "grey78"

  # Horizontal zero line always drawn: y = OLS residual is always centred.
  # Vertical reference: zero when centred; within-group mean of x when not.
  # The group mean falls in the centre of each column of panels and gives
  # the same above/below average split as zero does in the centred case.
  vline_layer <- if (centre) {
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30")
  } else {
    x_means <- aggregate(x ~ Group, data = df, FUN = mean)
    ggplot2::geom_vline(
      data     = x_means,
      ggplot2::aes(xintercept = x),
      linetype = "dotted", linewidth = 0.7, colour = "grey30"
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                        linewidth  = 0.7, colour = "grey30") +
    vline_layer +
    ggplot2::geom_point(size = 1.6, alpha = 0.7, colour = base_col, ...) +
    .freg_highlight_layers(df, hl) +
    ggplot2::facet_grid(pair_label ~ Group, scales = "free") +
    ggplot2::labs(
      x       = if (centre) "Efficiency (centred BLUE)" else
                             "Efficiency (unconditional BLUE)",
      y       = "Response index (OLS residual)",
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


# ---- Public function -----------------------------------------------------

#' Plot Output from fixedRegress()
#'
#' @description
#' Produces one of two ggplot2 visualisations from a [fixedRegress()] result
#' object, analogous to [plot_randomRegress()] for fixed-effects (BLUE) output.
#' Both plot types return a `ggplot` object that can be further customised
#' with the standard `+` operator.
#'
#' @details
#' The two `type` options are:
#' \describe{
#'   \item{`"regress"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     group (columns).  Each panel plots the raw conditioned-treatment BLUEs
#'     (y) against the conditioning-treatment BLUEs (x) for one group × one
#'     treatment pair.  A dotted OLS regression line is drawn and the slope
#'     \eqn{\hat{\beta}} is annotated in the top-left corner of each panel.
#'     Unlike the random-effects analogue the line need not pass through the
#'     origin because OLS includes an intercept.}
#'   \item{`"quadrant"`}{Grid of scatter plots faceted by BLUP pair (rows) and
#'     group (columns).  Each panel plots the response index (OLS residual,
#'     y) against the unconditional treatment BLUE (efficiency, x).  Dotted
#'     zero reference lines divide each panel into four quadrants.}
#' }
#'
#' **Variety highlighting** — see [plot_randomRegress()] for the full
#' description of the polar-angle algorithm used to select the default six
#' genotypes.  The same three archetypes (efficiency, balanced,
#' responsiveness) are chosen from the top-right and bottom-left quadrants of
#' the efficiency × response-index space, averaged across all groups.
#'
#' The output data frames returned when `return_data = TRUE` use standardised
#' column names `"Group"` (the grouping variable) and `"Genotype"` (the
#' regression unit) regardless of the original column names in
#' `res$blues`.
#'
#' @param res         A list returned by [fixedRegress()].
#' @param type        Character string selecting the plot type. One of
#'   `"regress"` (default) or `"quadrant"`.
#' @param treatments  Character vector restricting which conditioned treatments
#'   are included.  `NULL` (default) includes all conditioned treatments.
#' @param centre      Logical.  If `TRUE` (default), within-group means are
#'   subtracted from the x-axis (unconditional treatment BLUE) and, for
#'   `type = "regress"`, also from the y-axis.  This removes the treatment
#'   fixed-effect mean so that the zero reference lines fall within the data
#'   and the quadrant concept is meaningful.  Set to `FALSE` to display raw
#'   BLUEs including the treatment mean.
#' @param highlight   Controls genotype annotation. One of:
#'   \describe{
#'     \item{`"default"`}{Automatically selects 6 genotypes using the
#'       polar-angle algorithm (3 top-right, 3 bottom-left).}
#'     \item{character vector}{Highlight exactly the named genotypes.}
#'     \item{`NULL`}{No highlighting.}
#'   }
#'   Default is `"default"`.
#' @param theme       A complete ggplot2 theme object. Default
#'   [ggplot2::theme_bw()].
#' @param return_data Logical. If `TRUE` returns the tidy data frame rather
#'   than the plot. Default `FALSE`.
#' @param ...         Additional arguments passed to the background
#'   `geom_point()` call (e.g. `size`, `alpha`).
#'
#' @return A `ggplot` object (when `return_data = FALSE`) or a `data.frame`
#'   (when `return_data = TRUE`).
#'
#' @seealso [fixedRegress()], [plot_randomRegress()], [ggplot2::ggplot()]
#'
#' @examples
#' \dontrun{
#' res <- fixedRegress(model, term = "Treatment:Site:Genotype",
#'                     by = "Site", levs = c("T0", "T1", "T2"))
#'
#' # Regression plot with default 6 highlighted genotypes
#' plot_fixedRegress(res)
#'
#' # Quadrant plot
#' plot_fixedRegress(res, type = "quadrant")
#'
#' # User-specified highlights
#' plot_fixedRegress(res, type = "quadrant",
#'                   highlight = c("Gen01", "Gen15"))
#'
#' # Suppress highlighting
#' plot_fixedRegress(res, type = "regress", highlight = NULL)
#'
#' # Retrieve tidy data frame
#' df <- plot_fixedRegress(res, type = "quadrant", return_data = TRUE)
#' }
#'
#' @export
plot_fixedRegress <- function(res,
                              type        = c("regress", "quadrant"),
                              treatments  = NULL,
                              highlight   = "default",
                              centre      = TRUE,
                              theme       = ggplot2::theme_bw(),
                              return_data = FALSE,
                              ...) {

  # ---- Validate ----------------------------------------------------------
  type <- match.arg(type)

  required <- c("blues", "beta", "sigmat", "cond_list", "type")
  missing  <- setdiff(required, names(res))
  if (length(missing))
    stop("'res' is missing expected element(s): ",
         paste(missing, collapse = ", "),
         ". Was it produced by fixedRegress()?")

  if (!is.null(treatments) && !is.character(treatments))
    stop("'treatments' must be a character vector or NULL.")

  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object, e.g. ggplot2::theme_bw().")

  if (!is.logical(centre) || length(centre) != 1L)
    stop("'centre' must be a single logical value (TRUE or FALSE).")

  # ---- Build tidy data ---------------------------------------------------
  df <- switch(type,
    regress  = .freg_regress_data( res, treatments, centre),
    quadrant = .freg_quadrant_data(res, treatments, centre)
  )

  if (return_data) return(df)

  # ---- Resolve highlights ------------------------------------------------
  # Highlights are ALWAYS derived from centred quadrant data regardless of
  # the display 'centre' setting.  Without centring, x = absolute BLUE
  # (~4500) is always positive so the bottom-left quadrant is empty and no
  # blue archetypes can be selected.
  hl <- NULL
  if (identical(highlight, "default")) {
    qdata <- .freg_quadrant_data(res, treatments, centre = TRUE)
    hl    <- .freg_default_highlights(qdata)
  } else if (is.character(highlight) && length(highlight) > 0L) {
    all_genos <- unique(as.character(df$Genotype))
    bad       <- setdiff(highlight, all_genos)
    if (length(bad))
      warning("'highlight' genotypes not found in data: ",
              paste(bad, collapse = ", "))
    valid <- intersect(highlight, all_genos)
    if (length(valid) > 0L)
      hl <- data.frame(Genotype = valid, group = "custom",
                       stringsAsFactors = FALSE)
  }

  # ---- Build plot --------------------------------------------------------
  switch(type,
    regress  = .freg_plot_regress( df, hl, centre, theme, ...),
    quadrant = .freg_plot_quadrant(df, hl, centre, theme, ...)
  )
}
