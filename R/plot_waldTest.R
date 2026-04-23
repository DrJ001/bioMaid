# ============================================================
#  plot_waldTest.R
#  Forest plot visualisation for waldTest() output.
#
#  Only Contrast results can be displayed — a forest plot
#  requires an estimate and standard error, which joint zero
#  tests do not provide.  Supplying a result that contains
#  only Zero tests raises an error; supplying one that also
#  contains Contrasts raises a warning and the Zero results
#  are silently dropped.
#
#  Design:
#   - Contrast rows : filled circles with horizontal CI bars.
#   - Colour        : -log10(p) gradient — grey (non-sig) ->
#                     gold -> red (strongly significant).
#   - p-value label : printed beside each row.
#   - by-group      : facet_wrap(ncol=1, free_y) OR single panel
#                     with dotted group separators.
#
#  Pattern mirrors plot_fixedRegress / plot_randomRegress:
#    public function -> .wt_build_data() -> .wt_build_plot()
# ============================================================


# ---- Private helpers -------------------------------------------------------

## Detect the grouping column added by waldTest() when `by` was supplied.
## Returns the column name as a character string, or NULL.
#' @noRd
.wt_detect_by <- function(res) {
  known_con <- c("Comparison", "Estimate", "Std.Error",
                 "Wald.Statistic", "F.Statistic", "df", "P.Value")
  by_col <- NULL
  if (!is.null(res$Contrasts)) {
    extra <- setdiff(names(res$Contrasts), known_con)
    if (length(extra) == 1L) by_col <- extra
  }
  by_col
}

## Rescale a numeric vector to [0, 1].  Used to build gradientn breakpoints
## without requiring the 'scales' package.
#' @noRd
.wt_rescale01 <- function(x) {
  rng <- range(x, finite = TRUE)
  if (diff(rng) < .Machine$double.eps) return(rep(0.5, length(x)))
  (x - rng[1L]) / (rng[2L] - rng[1L])
}

## Build the tidy data frame used for plotting (Contrasts only).
#' @noRd
.wt_build_data <- function(res, ci_level, alpha, by_col) {

  z_mult   <- stats::qnorm(1 - (1 - ci_level) / 2)
  stat_col <- if (res$test == "F") "F.Statistic" else "Wald.Statistic"

  d     <- res$Contrasts
  ci_lo <- d$Estimate - z_mult * d$Std.Error
  ci_hi <- d$Estimate + z_mult * d$Std.Error

  out <- data.frame(
    label       = as.character(d$Comparison),
    group       = if (!is.null(by_col)) as.character(d[[by_col]]) else "All",
    Estimate    = d$Estimate,
    Std.Error   = d$Std.Error,
    CI_lower    = ci_lo,
    CI_upper    = ci_hi,
    x_label     = ci_hi,   # p-value label placed at right CI arm
    statistic   = d[[stat_col]],
    df_test     = d$df,
    P.Value     = d$P.Value,
    neg_log10_p = pmin(-log10(pmax(d$P.Value, .Machine$double.eps)), 50),
    significant = d$P.Value < alpha,
    stringsAsFactors = FALSE
  )

  # p-value display label
  out$p_label <- ifelse(
    out$P.Value < 0.001,
    "p<0.001",
    paste0("p=", formatC(out$P.Value, digits = 3L, format = "f"))
  )

  # Unique row ID prevents label collision when the same contrast label
  # appears in multiple groups (e.g. "N1 vs N0" in both Site1 and Site2).
  out$row_id <- paste(out$group, out$label, sep = "\u2060")

  # y-axis factor: reverse so the first row in `out` appears at the top.
  out$y_factor <- factor(out$row_id, levels = rev(unique(out$row_id)))

  out
}

## Build the ggplot forest plot object.
#' @noRd
.wt_build_plot <- function(data, res, facet, alpha, by_col, ci_level,
                            theme, ...) {

  # ---- Colour ramp --------------------------------------------------------
  # Non-significant: grey.  Significant: gold -> orange -> red -> dark red.
  # The break point is at -log10(alpha).
  nlp_thresh <- -log10(alpha + .Machine$double.eps)
  nlp_max    <- max(data$neg_log10_p, na.rm = TRUE)
  nlp_max    <- max(nlp_max, nlp_thresh * 1.5, 2.0)  # ensure visible range

  if (nlp_max <= nlp_thresh + 0.05) {
    # All results non-significant — grey gradient only
    colour_vals <- c(0, 1)
    colour_cols <- c("grey85", "grey42")
  } else {
    # Two-zone ramp: grey up to threshold, warm gradient above
    bp <- c(0,
            nlp_thresh * 0.97,
            nlp_thresh,
            nlp_thresh + (nlp_max - nlp_thresh) * 0.33,
            nlp_thresh + (nlp_max - nlp_thresh) * 0.67,
            nlp_max)
    colour_vals <- .wt_rescale01(bp)
    colour_cols <- c("grey82", "grey58", "#FDB863", "#E66101",
                     "#B2182B", "#67001F")
  }

  # ---- Labels & flags -----------------------------------------------------
  ci_pct   <- paste0(round(ci_level * 100L), "%")
  stat_lab <- if (res$test == "F") "F" else "Wald"
  adj_lab  <- if (!is.null(res$adjust) && res$adjust != "none")
    paste0("  |  p-adj: ", res$adjust) else ""

  has_grps <- !is.null(by_col) && length(unique(data$group)) > 1L

  # Map row_id -> readable label for scale_y_discrete
  y_labs <- stats::setNames(data$label, data$row_id)

  # ---- Base plot ----------------------------------------------------------
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = Estimate, y = y_factor, colour = neg_log10_p)
  ) +

    # Vertical reference line at x = 0
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      colour     = "grey30",
      linewidth  = 0.45
    ) +

    # CI bars
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
      height    = 0.22,
      linewidth = 0.9
    ) +

    # Point estimates
    ggplot2::geom_point(size = 3) +

    # p-value labels placed just right of the upper CI arm
    ggplot2::geom_text(
      ggplot2::aes(x = x_label, label = p_label),
      hjust       = -0.12,
      vjust       = 0.4,
      size        = 2.7,
      show.legend = FALSE
    ) +

    # ---- Scales -----------------------------------------------------------
    ggplot2::scale_y_discrete(labels = y_labs) +

    ggplot2::scale_colour_gradientn(
      colours  = colour_cols,
      values   = colour_vals,
      limits   = c(0, nlp_max),
      name     = expression(-log[10](italic(p))),
      guide    = ggplot2::guide_colourbar(
        title.position = "top",
        barwidth       = ggplot2::unit(8,    "lines"),
        barheight      = ggplot2::unit(0.65, "lines")
      )
    ) +

    # Expand x to the right to accommodate p-value text labels
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.28))
    ) +

    ggplot2::labs(
      x       = paste0("Estimate  [", ci_pct, " CI]"),
      y       = NULL,
      caption = paste0(stat_lab, "-test", adj_lab)
    ) +

    theme +

    ggplot2::theme(
      legend.position    = "top",
      legend.box         = "horizontal",
      panel.grid.major.y = ggplot2::element_line(colour    = "grey93",
                                                  linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_line(colour    = "grey88",
                                                  linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 9),
      strip.background   = ggplot2::element_rect(fill   = "grey92",
                                                  colour = NA),
      strip.text         = ggplot2::element_text(face = "bold", size = 9),
      plot.caption       = ggplot2::element_text(size   = 7.5,
                                                  colour = "grey48",
                                                  hjust  = 0)
    )

  # ---- Faceting -----------------------------------------------------------
  if (facet && has_grps) {
    # One panel per group; free_y so each panel shows only its own rows
    p <- p + ggplot2::facet_wrap(
      ~ group,
      ncol           = 1L,
      scales         = "free_y",
      strip.position = "top"
    )

  } else if (!facet && has_grps) {
    # Single panel: add dotted horizontal separator lines between groups.
    # Within the y factor, level 1 is at the bottom and level n at the top;
    # as.numeric(y_factor) gives each row's level number.
    # The separator sits 0.5 units below the lowest y-position of each group
    # (except the bottom group which needs no separator below it).
    grp_min_y  <- tapply(as.numeric(data$y_factor), data$group, min)
    bottom_grp <- names(which.min(grp_min_y))
    sep_y      <- grp_min_y[names(grp_min_y) != bottom_grp] - 0.5

    if (length(sep_y)) {
      p <- p + ggplot2::geom_hline(
        yintercept = unname(sep_y),
        linetype   = "dotted",
        colour     = "grey52",
        linewidth  = 0.4
      )
    }
  }

  p
}


# ---- Main exported function ------------------------------------------------

#' Forest Plot for Wald Test Results
#'
#' @description
#' Produces a forest plot from the `$Contrasts` component of a [waldTest()]
#' result.  A forest plot requires an estimate and standard error for each row,
#' so joint zero tests (which have neither) cannot be displayed.  If `res`
#' contains only zero-test results an error is raised; if it contains both
#' contrast and zero-test results a warning is issued and the zero results are
#' dropped.
#'
#' Each contrast is shown as a filled circle with horizontal confidence interval
#' bars.  Points are coloured by \eqn{-\log_{10}(p)}: non-significant results
#' (p \eqn{\geq} `alpha`) appear in grey; significant results follow a warm
#' gradient from gold through to dark red as evidence strengthens.  The raw
#' p-value is printed beside each row.  A vertical dashed line marks
#' \eqn{x = 0}.
#'
#' @param res         List returned by [waldTest()].
#' @param facet       Logical.  When [waldTest()] was called with a `by`
#'   argument, `TRUE` (default) produces one facet panel per group with free
#'   y-scales; `FALSE` places all groups on a single panel separated by dotted
#'   horizontal lines.
#' @param ci_level    Numeric in (0, 1).  Confidence level for the CI arms.
#'   Default `0.95`.
#' @param alpha       Numeric significance threshold.  Controls the colour-scale
#'   break between non-significant (grey) and significant (warm gradient)
#'   results, and the `significant` column in the returned data frame.
#'   Default `0.05`.
#' @param theme       A ggplot2 theme object used as the base theme.
#'   Default [ggplot2::theme_bw()].
#' @param return_data Logical.  If `TRUE`, returns the tidy data frame used to
#'   build the plot rather than the ggplot object.  Useful for bespoke
#'   customisation.  Default `FALSE`.
#' @param ...         Reserved for future use.
#'
#' @return A [ggplot2::ggplot] object returned invisibly (so it can be extended
#'   with `+`), or a data frame when `return_data = TRUE`.
#'
#' @examples
#' \dontrun{
#' ## Pairwise contrasts, no by-group
#' res <- waldTest(pred,
#'   cc = list(list(coef = c("N0","N1","N2"), type = "con", comp = "pairwise")))
#' plot_waldTest(res)
#'
#' ## By-group — one panel per site (default)
#' res2 <- waldTest(pred_site,
#'   cc = list(list(coef = c("N0","N1","N2"), type = "con", comp = "pairwise")),
#'   by = "Site")
#' plot_waldTest(res2, facet = TRUE)
#'
#' ## By-group — all sites on one panel with separator lines
#' plot_waldTest(res2, facet = FALSE)
#'
#' ## 99% CI
#' plot_waldTest(res, ci_level = 0.99)
#'
#' ## Stricter significance threshold
#' plot_waldTest(res, alpha = 0.01)
#'
#' ## Return tidy data for bespoke use
#' df <- plot_waldTest(res, return_data = TRUE)
#' }
#'
#' @seealso [waldTest()]
#' @export
plot_waldTest <- function(res,
                          facet       = TRUE,
                          ci_level    = 0.95,
                          alpha       = 0.05,
                          theme       = ggplot2::theme_bw(),
                          return_data = FALSE,
                          ...) {

  # ---- Validate res -------------------------------------------------------
  if (!is.list(res) || !all(c("test", "adjust") %in% names(res)))
    stop("'res' must be the list returned by waldTest().")

  # Zero-only: no point estimates or SEs — cannot draw a forest plot.
  if (is.null(res$Contrasts) && !is.null(res$Zero))
    stop("'res' contains only joint zero test results. ",
         "A forest plot requires contrast estimates and standard errors. ",
         "Run waldTest() with 'con' type contrasts to use plot_waldTest().")

  # No results at all.
  if (is.null(res$Contrasts))
    stop("'res' contains no Contrast results to plot.")

  # Contrasts present but Zero results also supplied — warn and proceed.
  if (!is.null(res$Zero))
    warning("Joint zero test results in 'res$Zero' cannot be displayed in a ",
            "forest plot (no estimates or standard errors) and have been ",
            "ignored. Only 'res$Contrasts' will be plotted.",
            call. = FALSE)

  # ---- Validate arguments -------------------------------------------------
  if (!is.logical(facet) || length(facet) != 1L)
    stop("'facet' must be TRUE or FALSE.")
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      ci_level <= 0 || ci_level >= 1)
    stop("'ci_level' must be a single number strictly between 0 and 1.")
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a single number strictly between 0 and 1.")
  if (!inherits(theme, "theme"))
    stop("'theme' must be a ggplot2 theme object.")
  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("'return_data' must be TRUE or FALSE.")

  # ---- Detect by-group column ---------------------------------------------
  by_col <- .wt_detect_by(res)

  # ---- Build tidy data frame ----------------------------------------------
  data <- .wt_build_data(res,
                         ci_level = ci_level,
                         alpha    = alpha,
                         by_col   = by_col)

  if (return_data) return(data)

  # ---- Build and return ggplot --------------------------------------------
  p <- .wt_build_plot(data,
                      res      = res,
                      facet    = facet,
                      alpha    = alpha,
                      by_col   = by_col,
                      ci_level = ci_level,
                      theme    = theme,
                      ...)
  invisible(p)
}
