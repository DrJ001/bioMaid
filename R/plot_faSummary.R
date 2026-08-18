# Non-standard-evaluation column names used in ggplot2 aesthetics
utils::globalVariables(c(
  "x_num", "ymin", "ymax", "source", "env_row", "env_col", "corr",
  "xend", "yend", "radius", "within_tol", "lab", "pair", "loading",
  "blup_adj", "slope", "facet_lab", "circ_x", "circ_y"
))

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Resolve which term of a faSummary object to plot
#' @noRd
.pfs_term <- function(res, term, need_loads = TRUE) {
  if (is.null(term)) {
    # Default to the first term that carries loadings
    cand <- if (need_loads) {
      Filter(function(nm) !is.null(res$gammas[[nm]]$loads), res$terms)
    } else {
      res$terms
    }
    if (length(cand) == 0L)
      stop("plot_faSummary(): no suitable FA term found in 'res'.", call. = FALSE)
    term <- cand[1L]
  }

  if (!term %in% res$terms)
    stop("plot_faSummary(): term '", term, "' not found in 'res'.\n  Available: ",
         paste(res$terms, collapse = ", "), call. = FALSE)

  g <- res$gammas[[term]]
  if (need_loads && is.null(g$loads))
    stop("plot_faSummary(): term '", term, "' has no loadings. ",
         "The combined vm() + ide() structure is not itself a factor analytic ",
         "model, so only type = \"heatmap\" is available for it.", call. = FALSE)

  list(term = term, g = g, b = res$blups[[term]])
}

#' Environment ordering for the heatmap and VAF plots
#' @noRd
.pfs_env_order <- function(g, order) {
  envs <- rownames(g$Cmat)
  switch(order,
    asis    = envs,
    loading = {
      if (is.null(g$loads)) envs else envs[order(-g$loads[, 1L])]
    },
    cluster = {
      hc <- stats::hclust(stats::as.dist(1 - g$Cmat), method = "average")
      envs[hc$order]
    }
  )
}

#' Non-overlapping text layer, falling back to geom_text without ggrepel
#' @noRd
.pfs_text_layer <- function(mapping, data, size = 2.6, ...) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(mapping = mapping, data = data, size = size,
                             max.overlaps = Inf, ...)
  } else {
    ggplot2::geom_text(mapping = mapping, data = data, size = size, ...)
  }
}

#' Score matrix (varieties x k) from the scores data frame
#' @noRd
.pfs_score_matrix <- function(g, b) {
  if (is.null(b) || is.null(b$scores))
    stop("plot_faSummary(): factor scores not available. ",
         "Re-run faSummary() with blups = TRUE.", call. = FALSE)
  sc  <- b$scores
  mat <- matrix(sc$blupr, ncol = g$k)
  rownames(mat) <- unique(as.character(sc[[g$inner]]))
  colnames(mat) <- paste0("fac_", seq_len(g$k))
  mat
}

#' Default variety selection: the extremes of each factor
#' @noRd
.pfs_default_varieties <- function(score_mat, n) {
  nms <- rownames(score_mat)
  picks <- unlist(lapply(seq_len(ncol(score_mat)), function(r) {
    c(nms[which.min(score_mat[, r])], nms[which.max(score_mat[, r])])
  }))
  picks <- unique(picks)
  if (length(picks) > n) picks <- picks[seq_len(n)]
  picks
}

# ---------------------------------------------------------------------------
# "VAF"
# ---------------------------------------------------------------------------

#' @noRd
.pfs_plot_VAF <- function(g, order, theme) {
  vaf_env  <- g$vaf_env
  vaf_summ <- g$vaf_summary
  k        <- g$k
  sterm    <- g$outer
  factor_cols <- paste0("Factor", seq_len(k))

  env_ord <- .pfs_env_order(g, order)
  src_levels <- c(paste0("Factor ", seq_len(k)), "Specific")

  # Pre-compute cumulative positions so each segment is placed exactly where
  # it belongs, regardless of near-zero proportions.  Factor 1 sits at the
  # bottom (ymin = 0) and Specific at the top (ymax = 1).  This avoids
  # position_stack entirely, which is prone to artefacts when a near-zero
  # segment is sandwiched between larger ones.
  long_rows <- vector("list", nrow(vaf_env) * (k + 1L))
  idx <- 1L
  for (i in seq_len(nrow(vaf_env))) {
    props <- c(unlist(vaf_env[i, factor_cols, drop = TRUE]),
               Specific = vaf_env$Specific[i])
    cumul <- c(0, cumsum(props))
    for (s in seq_along(src_levels)) {
      long_rows[[idx]] <- data.frame(
        env_fac = as.character(vaf_env[[sterm]][i]),
        source  = src_levels[s],
        ymin    = cumul[s],
        ymax    = cumul[s + 1L],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  long_df <- do.call(rbind, long_rows)
  long_df$env_fac <- factor(long_df$env_fac, levels = env_ord)
  long_df$source  <- factor(long_df$source,  levels = src_levels)
  long_df$x_num   <- as.integer(long_df$env_fac)

  pal <- c(
    if (k == 1L) "#2166AC"
    else grDevices::colorRampPalette(c("#2166AC", "#41B6C4"))(k),
    "grey75"
  )
  names(pal) <- src_levels

  spec_overall      <- vaf_summ$pct_var[vaf_summ$factor == "Specific"]
  overall_explained <- 1 - spec_overall
  bar_w <- 0.41

  plt <- ggplot2::ggplot(long_df) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = x_num - bar_w, xmax = x_num + bar_w,
                   ymin = ymin, ymax = ymax, fill = source),
      colour = NA
    ) +
    ggplot2::geom_hline(yintercept = overall_explained, linetype = "dashed",
                        colour = "grey30", linewidth = 0.5) +
    ggplot2::annotate("text", x = length(env_ord) + 0.4,
                      y = overall_explained,
                      label = sprintf("Mean\n%.0f%%", 100 * overall_explained),
                      hjust = 0, vjust = 0.5, size = 2.8, colour = "grey30") +
    ggplot2::scale_fill_manual(values = pal, name = "Source",
                               guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_x_continuous(breaks = seq_along(env_ord), labels = env_ord,
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = sterm, y = "Proportion of genetic variance",
      title = "Variance Accounted For (VAF) by FA factors per environment",
      subtitle = sprintf(
        "Overall: %.1f%% explained by FA factors  |  %.1f%% specific (unexplained)",
        100 * overall_explained, 100 * spec_overall)
    ) +
    theme +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      plot.subtitle   = ggplot2::element_text(size = 8, colour = "grey40"),
      legend.key.size = ggplot2::unit(0.45, "cm"),
      plot.margin     = ggplot2::margin(5, 40, 5, 5)
    )

  list(plot = plt, data = list(env = vaf_env, summary = vaf_summ))
}

# ---------------------------------------------------------------------------
# "heatmap"
# ---------------------------------------------------------------------------

#' @noRd
.pfs_plot_heatmap <- function(g, order, theme) {
  env_ord <- .pfs_env_order(g, order)
  cmat    <- g$Cmat[env_ord, env_ord, drop = FALSE]

  df <- data.frame(
    env_row = factor(rep(env_ord, times = length(env_ord)), levels = env_ord),
    env_col = factor(rep(env_ord, each  = length(env_ord)),
                     levels = rev(env_ord)),
    corr    = as.vector(cmat),
    stringsAsFactors = FALSE
  )

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = env_row, y = env_col,
                                          fill = corr)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), name = "Genetic\ncorrelation"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = g$outer, y = g$outer,
      title = "Genetic correlation between environments",
      subtitle = paste0("Ordering: ", order)
    ) +
    theme +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
      axis.text.y   = ggplot2::element_text(size = 7),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40")
    )

  list(plot = plt, data = df)
}

# ---------------------------------------------------------------------------
# "loadings"
# ---------------------------------------------------------------------------

#' @noRd
.pfs_plot_loadings <- function(g, tol, theme) {
  lam  <- g$loads_cor
  k    <- g$k
  envs <- rownames(lam)

  pairs <- utils::combn(seq_len(k), 2L)
  rows  <- vector("list", ncol(pairs))
  for (p in seq_len(ncol(pairs))) {
    i <- pairs[1L, p]
    j <- pairs[2L, p]
    rows[[p]] <- data.frame(
      pair   = sprintf("Factor %d vs Factor %d", i, j),
      lab    = envs,
      xend   = lam[, i],
      yend   = lam[, j],
      radius = sqrt(lam[, i]^2 + lam[, j]^2),
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  df$within_tol <- df$radius >= tol
  rownames(df)  <- NULL

  # Unit circle, repeated once per panel
  ang  <- seq(0, 2 * pi, length.out = 200L)
  circ <- do.call(rbind, lapply(unique(df$pair), function(pp) {
    data.frame(pair = pp, circ_x = cos(ang), circ_y = sin(ang),
               stringsAsFactors = FALSE)
  }))

  plt <- ggplot2::ggplot() +
    ggplot2::geom_path(data = circ,
                       ggplot2::aes(x = circ_x, y = circ_y),
                       colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_segment(
      data = df,
      ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend,
                   colour = within_tol, linetype = within_tol),
      linewidth = 0.5
    ) +
    .pfs_text_layer(
      mapping = ggplot2::aes(x = xend, y = yend, label = lab,
                             colour = within_tol),
      data = df, size = 2.4, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
      labels = c(`TRUE` = paste0("radius >= ", tol),
                 `FALSE` = paste0("radius < ", tol)),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = "none"
    ) +
    ggplot2::coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
    ggplot2::facet_wrap(~ pair) +
    ggplot2::labs(
      x = "Loading (first factor of pair)",
      y = "Loading (second factor of pair)",
      title = "Correlation-scaled FA loadings",
      subtitle = paste0(
        "Environments near the unit circle are well explained by the pair; ",
        "angle between vectors approximates genetic correlation")
    ) +
    theme +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      strip.text    = ggplot2::element_text(size = 8)
    )

  list(plot = plt, data = df)
}

# ---------------------------------------------------------------------------
# "regress" and "added"
# ---------------------------------------------------------------------------

#' @noRd
.pfs_plot_regress <- function(g, b, varieties, added, theme) {
  if (is.null(b) || is.null(b$blups))
    stop("plot_faSummary(): BLUPs not available. ",
         "Re-run faSummary() with blups = TRUE.", call. = FALSE)

  lam       <- g$loads
  k         <- g$k
  score_mat <- .pfs_score_matrix(g, b)

  missing_v <- setdiff(varieties, rownames(score_mat))
  if (length(missing_v) > 0L)
    stop("plot_faSummary(): variety not found: ",
         paste(missing_v, collapse = ", "), call. = FALSE)

  bl <- b$blups
  bl <- bl[bl[[g$inner]] %in% varieties, , drop = FALSE]

  rows <- vector("list", length(varieties) * k)
  idx  <- 1L
  for (v in varieties) {
    bv  <- bl[bl[[g$inner]] == v, , drop = FALSE]
    yv  <- bv$blup[match(rownames(lam), bv[[g$outer]])]
    for (r in seq_len(k)) {
      # Added-variable adjustment: remove the contribution of every other
      # factor so the plotted slope is the factor-r score.
      adj <- if (added && k > 1L) {
        drop(lam[, -r, drop = FALSE] %*% score_mat[v, -r])
      } else {
        rep(0, nrow(lam))
      }
      rows[[idx]] <- data.frame(
        variety   = v,
        facet_lab = paste0("Factor ", r),
        env       = rownames(lam),
        loading   = lam[, r],
        blup_adj  = yv - adj,
        slope     = score_mat[v, r],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  df <- do.call(rbind, rows)
  df$variety <- factor(df$variety, levels = varieties)
  rownames(df) <- NULL

  slopes <- unique(df[, c("variety", "facet_lab", "slope")])

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = loading, y = blup_adj)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_abline(
      data = slopes,
      ggplot2::aes(slope = slope, intercept = 0),
      colour = "#2166AC", linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 1.1, alpha = 0.8) +
    ggplot2::facet_grid(variety ~ facet_lab) +
    ggplot2::labs(
      x = "Environment loading",
      y = if (added) "Adjusted BLUP" else "BLUP",
      title = if (added)
        "Added variable plot: BLUP vs loading, other factors removed"
      else
        "Regression of BLUPs on environment loadings",
      subtitle = "Line slope is the variety's factor score"
    ) +
    theme +
    ggplot2::theme(
      plot.subtitle    = ggplot2::element_text(size = 8, colour = "grey40"),
      strip.text       = ggplot2::element_text(size = 7),
      panel.spacing    = ggplot2::unit(0.3, "lines")
    )

  list(plot = plt, data = df)
}

# ---------------------------------------------------------------------------
# plot_faSummary()
# ---------------------------------------------------------------------------

#' Visualise a factor analytic model summary
#'
#' Diagnostic plots for the factor analytic structure returned by
#' [faSummary()]: variance accounted for, the genetic correlation heatmap,
#' correlation-scaled loadings, and the regression of genotype BLUPs on
#' environment loadings.
#'
#' @details
#' # Plot types
#'
#' \describe{
#'   \item{`"VAF"`}{Stacked 100% bar chart of the Variance Accounted For in
#'     each environment: one segment per FA factor (sequential blues, from the
#'     bottom) plus the specific variance (grey, on top).  The dashed line
#'     marks the overall mean proportion explained.}
#'   \item{`"heatmap"`}{Genetic correlation between environments, on a
#'     diverging palette centred at zero.  Use `order` to reveal block
#'     structure.}
#'   \item{`"loadings"`}{Correlation-scaled loadings plotted as vectors from
#'     the origin, one panel per pair of factors, with the unit circle for
#'     reference.  Environments whose vector reaches beyond `tol` are well
#'     explained by that pair of factors and are drawn solid; the rest are
#'     dotted.  Requires \eqn{k \ge 2}.}
#'   \item{`"regress"`}{For each selected variety and factor, the genotype
#'     BLUPs plotted against the environment loadings, with a line whose slope
#'     is that variety's factor score.}
#'   \item{`"added"`}{As `"regress"`, but the BLUPs are first adjusted by
#'     removing the fitted contribution of every other factor, so the plotted
#'     slope aligns with the scatter.  This is the correct display when
#'     \eqn{k \ge 2}.}
#' }
#'
#' @param res         An object of class `"faSummary"` from [faSummary()].
#' @param type        Character; the visualisation to produce.  One of
#'   `"VAF"` (default), `"heatmap"`, `"loadings"`, `"regress"`, `"added"`.
#'   May be abbreviated.
#' @param term        Character; which FA term of `res` to plot.  The default,
#'   `NULL`, uses the first term that carries loadings.
#' @param order       Environment ordering for `"VAF"` and `"heatmap"`.  One of
#'   `"loading"` (default; by first-factor loading, descending), `"asis"` (the
#'   order of the factor levels) or `"cluster"` (average-linkage hierarchical
#'   clustering of \eqn{1 - C}).
#' @param varieties   Varieties to display in `"regress"` and `"added"`.
#'   Either `"default"` (the highest- and lowest-scoring variety on each
#'   factor) or a character vector of variety names.
#' @param n_varieties Positive integer; maximum number of varieties to select
#'   when `varieties = "default"`.  Default `6L`.
#' @param tol         Numeric in \[0, 1\]; radius threshold used to flag
#'   well-explained environments in the `"loadings"` plot.  Default `0.85`.
#' @param theme       A complete ggplot2 theme object.  Default
#'   [ggplot2::theme_bw()].
#' @param return_data Logical; if `TRUE`, return a list with elements `plot`
#'   and `data`.  Default `FALSE`.
#' @param ...         Reserved for future use; currently ignored.
#'
#' @return When `return_data = FALSE` (default) a `ggplot` object; otherwise a
#'   named list with elements `plot` and `data`.  For `type = "VAF"` the `data`
#'   element is itself a list with components `env` and `summary`.
#'
#' @seealso [faSummary()], [fastIC()], [plot_fastIC()]
#'
#' @examples
#' \dontrun{
#' fas <- faSummary(model)
#'
#' plot_faSummary(fas, type = "VAF")
#' plot_faSummary(fas, type = "heatmap", order = "cluster")
#' plot_faSummary(fas, type = "loadings")
#' plot_faSummary(fas, type = "added", varieties = c("Var01", "Var22"))
#'
#' out <- plot_faSummary(fas, type = "VAF", return_data = TRUE)
#' head(out$data$env)
#' }
#'
#' @export
plot_faSummary <- function(res,
                           type        = c("VAF", "heatmap", "loadings",
                                           "regress", "added"),
                           term        = NULL,
                           order       = c("loading", "asis", "cluster"),
                           varieties   = "default",
                           n_varieties = 6L,
                           tol         = 0.85,
                           theme       = ggplot2::theme_bw(),
                           return_data = FALSE,
                           ...) {

  # ---- Validate ----------------------------------------------------------
  if (!inherits(res, "faSummary"))
    stop("plot_faSummary(): 'res' must be an object of class 'faSummary'. ",
         "Was it produced by faSummary()?", call. = FALSE)

  type  <- match.arg(type)
  order <- match.arg(order)

  if (!inherits(theme, "theme"))
    stop("plot_faSummary(): 'theme' must be a ggplot2 theme object, ",
         "e.g. ggplot2::theme_bw().", call. = FALSE)

  if (!is.logical(return_data) || length(return_data) != 1L)
    stop("plot_faSummary(): 'return_data' must be a single logical value.",
         call. = FALSE)

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      tol < 0 || tol > 1)
    stop("plot_faSummary(): 'tol' must be a single number between 0 and 1.",
         call. = FALSE)

  n_varieties <- max(1L, as.integer(n_varieties[1L]))

  sel <- .pfs_term(res, term, need_loads = type != "heatmap")
  g   <- sel$g
  b   <- sel$b

  # ---- Type-specific preconditions ---------------------------------------
  if (type == "VAF" && is.null(g$vaf_env))
    stop("plot_faSummary(): VAF data not found for term '", sel$term, "'.",
         call. = FALSE)

  if (type == "loadings" && g$k < 2L)
    stop("plot_faSummary(): type = 'loadings' requires k >= 2 factors; ",
         "only ", g$k, " factor(s) found.", call. = FALSE)

  # ---- Resolve varieties for the regression plots ------------------------
  if (type %in% c("regress", "added")) {
    score_mat <- .pfs_score_matrix(g, b)
    varieties <- if (identical(varieties, "default")) {
      .pfs_default_varieties(score_mat, n_varieties)
    } else {
      as.character(varieties)
    }
    if (length(varieties) == 0L)
      stop("plot_faSummary(): no varieties selected.", call. = FALSE)
  }

  # ---- Dispatch ----------------------------------------------------------
  result <- switch(type,
    VAF      = .pfs_plot_VAF(g, order, theme),
    heatmap  = .pfs_plot_heatmap(g, order, theme),
    loadings = .pfs_plot_loadings(g, tol, theme),
    regress  = .pfs_plot_regress(g, b, varieties, added = FALSE, theme),
    added    = .pfs_plot_regress(g, b, varieties, added = TRUE,  theme)
  )

  if (return_data) result else result$plot
}
