# =============================================================================
# plot_accuracy.R  --  ggplot2 visualisations for accuracy() output
# =============================================================================
#
# Single-model types  (res2 not required):
#   "lollipop"  Group-level mean accuracy; lollipop + optional sd error bars.
#               Faceted by metric when both accuracy and gen.H2 are present.
#   "violin"    Per-variety accuracy distribution per group (requires by_variety
#               data). Faceted by metric when both are present.
#   "heatmap"   Variety x Group accuracy grid (requires by_variety data).
#               Faceted by metric when both are present.
#
# Comparison types  (res2 required):
#   "dumbbell"  Paired group-level dots connected by a coloured segment
#               (green = model2 improved, red = model2 declined).
#   "scatter"   Per-variety model1 (x) vs model2 (y) accuracy scatter with a
#               diagonal reference line (requires by_variety data in both).
#   "diff"      Signed accuracy gain (model2 - model1) lollipop per group.
# =============================================================================

utils::globalVariables(c(
  "group", "value", "sd_val", "metric_label", "model_label",
  "variety", "x1", "x2", "diff_val", "gain_dir",
  ".value_x", ".value_y",
  "present", "x_obs", "x_unobs", "pres_label"
))

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Plot accuracy results from \code{accuracy()}
#'
#' Provides six plot types for visualising BLUP accuracy computed by
#' \code{\link{accuracy}()}: a lollipop summary, violin distribution, heatmap,
#' and three two-model comparison plots (dumbbell, scatter, diff).
#'
#' @param res      Data frame returned by \code{accuracy()}. May be group-level
#'   or by-variety (\code{by_variety = TRUE}) depending on \code{type}.
#' @param res2     Optional second data frame from \code{accuracy()} for
#'   comparison plots (\code{"dumbbell"}, \code{"scatter"}, \code{"diff"}).
#' @param type     Plot type -- one of \code{"lollipop"}, \code{"violin"},
#'   \code{"heatmap"}, \code{"dumbbell"}, \code{"scatter"}, \code{"diff"}.
#' @param metric   Character vector: \code{"accuracy"}, \code{"gen.H2"}, or
#'   both (default). When both are selected the plot is faceted into two panels.
#' @param label1   Legend label for \code{res}  (comparison plots only).
#' @param label2   Legend label for \code{res2} (comparison plots only).
#' @param theme    A ggplot2 theme object applied to all plots.
#' @param return_data Logical. If \code{TRUE}, returns a list with elements
#'   \code{$plot} and \code{$data} instead of just the plot.
#' @param ...      Currently unused.
#'
#' @return A \code{ggplot} object, or a named list when \code{return_data = TRUE}.
#' @export
#'
#' @examples
#' \dontrun{
#' acc  <- accuracy(model_fa,   metric = c("accuracy","gen.H2"))
#' acc2 <- accuracy(model_diag, metric = c("accuracy","gen.H2"))
#'
#' plot_accuracy(acc,       type = "lollipop")
#' plot_accuracy(acc_bv,    type = "violin")
#' plot_accuracy(acc_bv,    type = "heatmap")
#' plot_accuracy(acc, acc2, type = "dumbbell", label1 = "FA", label2 = "Diag")
#' plot_accuracy(acc_bv, acc2_bv, type = "scatter")
#' plot_accuracy(acc, acc2, type = "diff", label1 = "FA", label2 = "Diag")
#' }
plot_accuracy <- function(res,
                           res2        = NULL,
                           type        = c("lollipop", "violin", "heatmap",
                                           "dumbbell", "scatter", "diff"),
                           metric      = c("accuracy", "gen.H2"),
                           label1      = "Model 1",
                           label2      = "Model 2",
                           theme       = ggplot2::theme_bw(),
                           return_data = FALSE,
                           ...) {

  type   <- match.arg(type)
  metric <- match.arg(metric, several.ok = TRUE)

  # ---- helpers ---------------------------------------------------------------
  is_bv  <- function(d) "variety" %in% names(d)

  # ---- type/data validation --------------------------------------------------
  single_types  <- c("lollipop", "violin", "heatmap")
  compare_types <- c("dumbbell", "scatter", "diff")
  bv_types      <- c("violin", "heatmap", "scatter")

  has_pres_spl <- "variety" %in% names(res) && "present" %in% names(res) &&
                  any(res$present) && any(!res$present)

  if (type %in% single_types && !is.null(res2))
    warning(sprintf("[plot_accuracy] type='%s' ignores res2.", type))
  if (type %in% compare_types && is.null(res2)) {
    if (type == "diff" && has_pres_spl) {
      # Allow: single-model present-gap diff
    } else {
      stop(sprintf("[plot_accuracy] type='%s' requires res2.", type))
    }
  }
  if (type %in% bv_types && !is_bv(res))
    stop(sprintf("[plot_accuracy] type='%s' requires by_variety=TRUE data in res.", type))
  if (type == "scatter" && !is.null(res2) && !is_bv(res2))
    stop("[plot_accuracy] type='scatter' requires by_variety=TRUE data in res2.")

  # ---- metric column validation ----------------------------------------------
  .check_metric <- function(d, m, lbl = "res") {
    acc_col <- if (is_bv(d)) "accuracy" else "mean_acc"
    if ("accuracy" %in% m && !acc_col %in% names(d))
      stop(sprintf(
        "[plot_accuracy] 'accuracy' not in %s -- run accuracy(..., metric='accuracy').", lbl))
    if ("gen.H2" %in% m && !"gen.H2" %in% names(d))
      stop(sprintf(
        "[plot_accuracy] 'gen.H2' not in %s -- run accuracy(..., metric='gen.H2').", lbl))
  }
  .check_metric(res, metric, "res")
  if (!is.null(res2)) .check_metric(res2, metric, "res2")

  # ---- dispatch --------------------------------------------------------------
  switch(type,
    lollipop = .pa_lollipop(res,           metric, theme, return_data),
    violin   = .pa_violin(  res,           metric, theme, return_data),
    heatmap  = .pa_heatmap( res,           metric, theme, return_data),
    dumbbell = .pa_dumbbell(res, res2, metric, label1, label2, theme, return_data),
    scatter  = .pa_scatter( res, res2, metric, label1, label2, theme, return_data),
    diff     = .pa_diff(    res, res2, metric, label1, label2, theme, return_data)
  )
}

# =============================================================================
# DATA PREP HELPERS
# =============================================================================

# Colour palette used throughout
.pa_metric_cols   <- c("Mrode Accuracy" = "#2166AC", "Gen. H\u00b2" = "#D6604D")
.pa_gain_cols     <- c("Improved" = "#1A9641", "Declined"  = "#D7191C")
.pa_present_cols  <- c("Observed" = "#2166AC", "Unobserved" = "#E07B39")

# Convert group-level accuracy() output to long format
.pa_long_grp <- function(res, metric) {
  out <- list()
  if ("accuracy" %in% metric && "mean_acc" %in% names(res)) {
    out[["accuracy"]] <- data.frame(
      group        = res$group,
      n_vars       = res$n_vars,
      value        = res$mean_acc,
      sd_val       = if ("sd_acc" %in% names(res)) res$sd_acc else NA_real_,
      metric_label = "Mrode Accuracy",
      stringsAsFactors = FALSE
    )
  }
  if ("gen.H2" %in% metric && "gen.H2" %in% names(res)) {
    out[["gen.H2"]] <- data.frame(
      group        = res$group,
      n_vars       = if ("n_vars" %in% names(res)) res$n_vars else NA_integer_,
      value        = res$gen.H2,
      sd_val       = NA_real_,
      metric_label = "Gen. H\u00b2",
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) stop("[plot_accuracy] No matching metric columns found in res.")
  do.call(rbind, c(out, make.row.names = FALSE))
}

# Convert by_variety accuracy() output to long format.
# The `present` column is carried through when it exists in res.
.pa_long_bv <- function(res, metric) {
  has_present <- "present" %in% names(res)
  out <- list()
  if ("accuracy" %in% metric && "accuracy" %in% names(res)) {
    df <- data.frame(
      group        = res$group,
      variety      = res$variety,
      value        = res$accuracy,
      metric_label = "Mrode Accuracy",
      stringsAsFactors = FALSE
    )
    if (has_present) df$present <- res$present
    out[["accuracy"]] <- df
  }
  if ("gen.H2" %in% metric && "gen.H2" %in% names(res)) {
    df <- data.frame(
      group        = res$group,
      variety      = res$variety,
      value        = res$gen.H2,
      metric_label = "Gen. H\u00b2",
      stringsAsFactors = FALSE
    )
    if (has_present) df$present <- res$present
    out[["gen.H2"]] <- df
  }
  if (!length(out)) stop("[plot_accuracy] No matching metric columns found in res.")
  do.call(rbind, c(out, make.row.names = FALSE))
}

# Compute group-level mean accuracy split by present status from by_variety data.
# Returns long data frame with columns: group, value, sd_val, metric_label,
# present (logical), pres_label (character).
.pa_present_grp <- function(res, metric) {
  out <- list()
  for (pres in c(TRUE, FALSE)) {
    sub <- res[res$present == pres, , drop = FALSE]
    if (nrow(sub) == 0L) next
    pl <- if (pres) "Observed" else "Unobserved"
    if ("accuracy" %in% metric && "accuracy" %in% names(sub)) {
      agg    <- aggregate(accuracy ~ group, data = sub, FUN = mean, na.rm = TRUE)
      agg_sd <- aggregate(accuracy ~ group, data = sub, FUN = sd,   na.rm = TRUE)
      out[[paste0("acc_", pres)]] <- data.frame(
        group        = agg$group,
        value        = agg$accuracy,
        sd_val       = agg_sd$accuracy,
        metric_label = "Mrode Accuracy",
        present      = pres,
        pres_label   = pl,
        stringsAsFactors = FALSE
      )
    }
    if ("gen.H2" %in% metric && "gen.H2" %in% names(sub)) {
      agg <- aggregate(gen.H2 ~ group, data = sub, FUN = function(x) x[1L])
      out[[paste0("h2_", pres)]] <- data.frame(
        group        = agg$group,
        value        = agg$gen.H2,
        sd_val       = NA_real_,
        metric_label = "Gen. H\u00b2",
        present      = pres,
        pres_label   = pl,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) stop("[plot_accuracy] No matching metric columns found in res.")
  do.call(rbind, c(out, make.row.names = FALSE))
}

# Collapse a by_variety accuracy() data frame to the group-level format that
# .pa_long_grp() expects (mean_acc, sd_acc, n_vars, gen.H2).
# Returns the input unchanged if it is already group-level (no 'variety' col).
.pa_bv_to_grp <- function(d) {
  if (!"variety" %in% names(d)) return(d)
  out <- list()
  if ("accuracy" %in% names(d)) {
    agg_m  <- aggregate(accuracy ~ group, data = d, FUN = mean, na.rm = TRUE)
    agg_s  <- aggregate(accuracy ~ group, data = d, FUN = sd,   na.rm = TRUE)
    agg_n  <- aggregate(accuracy ~ group, data = d,
                        FUN = function(x) sum(!is.na(x)))
    grp_df <- data.frame(
      group    = agg_m$group,
      mean_acc = agg_m$accuracy,
      sd_acc   = agg_s$accuracy,
      n_vars   = agg_n$accuracy,
      stringsAsFactors = FALSE
    )
    if ("gen.H2" %in% names(d)) {
      h2 <- aggregate(gen.H2 ~ group, data = d,
                      FUN = function(x) x[1L])
      grp_df <- merge(grp_df, h2, by = "group")
    }
    out <- grp_df
  }
  out
}

# Set group factor levels ordered by mean value (ascending -> best at top of y)
# Uses the first metric_label present if multiple exist.
.pa_order_groups <- function(d) {
  ref_metric <- unique(d$metric_label)[1L]
  ref        <- d[d$metric_label == ref_metric, ]
  ord        <- tapply(ref$value, ref$group, mean, na.rm = TRUE)
  d$group    <- factor(d$group, levels = names(sort(ord)))
  d
}

# Set variety factor levels ordered by mean value (descending -> best at top of y)
.pa_order_varieties <- function(d) {
  ref_metric <- unique(d$metric_label)[1L]
  ref        <- d[d$metric_label == ref_metric, ]
  ord        <- tapply(ref$value, ref$variety, mean, na.rm = TRUE)
  d$variety  <- factor(d$variety, levels = names(sort(ord, decreasing = TRUE)))
  d
}

# =============================================================================
# LOLLIPOP  (single model, group-level)
# =============================================================================

.pa_lollipop <- function(res, metric, theme, return_data) {

  # ---- Present-split mode: by_variety data with both observed and unobserved --
  is_bv        <- "variety" %in% names(res)
  has_pres_spl <- is_bv && "present" %in% names(res) &&
                  any(res$present) && any(!res$present)

  if (has_pres_spl) {
    d <- .pa_present_grp(res, metric)

    # Order groups by the observed (present=TRUE) mean accuracy, ascending
    ref_metric <- unique(d$metric_label)[1L]
    ref_obs    <- d[d$metric_label == ref_metric & d$present, ]
    grp_ord    <- names(sort(tapply(ref_obs$value, ref_obs$group, mean,
                                    na.rm = TRUE)))
    d$group    <- factor(d$group, levels = grp_ord)

    d_obs   <- d[d$present, ,  drop = FALSE]
    d_unobs <- d[!d$present, , drop = FALSE]

    # Join for the connecting segment
    seg <- merge(
      d_obs  [, c("group", "value", "metric_label")],
      d_unobs[, c("group", "value", "metric_label")],
      by = c("group", "metric_label"), suffixes = c("_obs", "_unobs")
    )
    names(seg)[names(seg) == "value_obs"]   <- "x_obs"
    names(seg)[names(seg) == "value_unobs"] <- "x_unobs"
    seg$group <- factor(seg$group, levels = levels(d$group))

    use_facet <- length(unique(d$metric_label)) > 1L

    p <- ggplot2::ggplot() +
      # Stem from 0 to observed value
      ggplot2::geom_segment(
        data = d_obs,
        ggplot2::aes(x = 0, xend = value, y = group, yend = group),
        colour = "grey65", linewidth = 0.7
      ) +
      # Connecting segment: observed to unobserved point
      ggplot2::geom_segment(
        data = seg,
        ggplot2::aes(x = x_obs, xend = x_unobs, y = group, yend = group),
        colour = "grey45", linewidth = 1.0, linetype = "dashed"
      ) +
      # Unobserved point (drawn first so observed sits on top)
      ggplot2::geom_point(
        data = d_unobs,
        ggplot2::aes(x = value, y = group, fill = pres_label),
        shape = 21, colour = "white", stroke = 0.4, size = 3.5
      ) +
      # Observed point
      ggplot2::geom_point(
        data = d_obs,
        ggplot2::aes(x = value, y = group, fill = pres_label),
        shape = 21, colour = "white", stroke = 0.4, size = 3.5
      ) +
      ggplot2::scale_fill_manual(
        name   = NULL,
        values = .pa_present_cols
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        expand = ggplot2::expansion(mult = c(0, 0.03))
      ) +
      ggplot2::labs(x = "Accuracy", y = NULL) +
      theme

    if (use_facet)
      p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2, scales = "free_x")

    if (return_data) return(list(plot = p, data = d))
    return(p)
  }

  # ---- Standard mode: group-level data ----------------------------------------
  d <- .pa_long_grp(res, metric)
  d <- .pa_order_groups(d)

  use_facet <- length(unique(d$metric_label)) > 1L

  p <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = group)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = value, yend = group),
      colour = "grey65", linewidth = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = metric_label),
      size = 3.5, show.legend = !use_facet
    )

  # Horizontal error bars for accuracy metric only (sd_val is NA for gen.H2)
  d_err <- d[!is.na(d$sd_val) & d$sd_val > 0, ]
  if (nrow(d_err) > 0L) {
    p <- p + ggplot2::geom_errorbar(
      data = d_err,
      ggplot2::aes(xmin = pmax(0, value - sd_val),
                   xmax = pmin(1, value + sd_val),
                   y    = group),
      orientation = "y",
      width    = 0.3,
      colour   = "grey40",
      linewidth = 0.5
    )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::scale_colour_manual(
      values = .pa_metric_cols,
      guide  = if (use_facet) "none" else ggplot2::guide_legend(title = NULL)
    ) +
    ggplot2::labs(x = "Accuracy", y = NULL, colour = NULL) +
    theme

  if (use_facet)
    p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2, scales = "free_x")

  if (return_data) return(list(plot = p, data = d))
  p
}

# =============================================================================
# VIOLIN  (single model, by_variety)
# =============================================================================

.pa_violin <- function(res, metric, theme, return_data) {

  d <- .pa_long_bv(res, metric)

  # Cullis H2 (gen.H2) is a group-level constant in by_variety data -- it
  # produces a degenerate flat display in a violin. Suppress it automatically
  # when Mrode Accuracy is also present; direct users to lollipop instead.
  H2_lbl     <- "Gen. H\u00b2"
  acc_lbl    <- "Mrode Accuracy"
  has_H2     <- H2_lbl  %in% unique(d$metric_label)
  has_acc    <- acc_lbl %in% unique(d$metric_label)

  if (has_H2 && has_acc) {
    message("[plot_accuracy] violin: Cullis H\u00b2 (gen.H2) is a group-level statistic ",
            "and is suppressed in violin mode. ",
            "Use type = \"lollipop\" to visualise it alongside Mrode Accuracy.")
    d <- d[d$metric_label == acc_lbl, ]
  }

  d <- .pa_order_groups(d)

  use_facet <- length(unique(d$metric_label)) > 1L

  # Detect if a metric's values are constant within every group (gen.H2 case
  # when only gen.H2 was requested and accuracy was absent)
  .is_const <- function(sub) {
    all(tapply(sub$value, sub$group,
               function(x) length(unique(round(x, 10)))) == 1L)
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = value))

  # When `present` is in the data, colour jitter points by observed status.
  has_present_col <- "present" %in% names(d)

  # Helper: build the jitter layer — fixed colour or coloured by present status
  .jitter_layer <- function(sub, col) {
    if (has_present_col && "present" %in% names(sub)) {
      sub$pres_label <- ifelse(sub$present, "Observed", "Unobserved")
      ggplot2::geom_jitter(
        data        = sub,
        ggplot2::aes(colour = pres_label),
        alpha = 0.75, size = 1.5,
        width = 0.10, height = 0,
        show.legend = TRUE
      )
    } else {
      ggplot2::geom_jitter(
        data   = sub,
        colour = col, alpha = 0.65, size = 1.5,
        width  = 0.10, height = 0
      )
    }
  }

  if (use_facet) {
    # Build per-metric geoms using data= subsets; facet handles layout
    for (ml in unique(d$metric_label)) {
      sub <- d[d$metric_label == ml, ]
      col <- .pa_metric_cols[ml]
      if (.is_const(sub)) {
        # Cullis H2 (gen.H2): one point per group (constant within group)
        grp_sum <- aggregate(value ~ group, data = sub, FUN = mean)
        p <- p +
          ggplot2::geom_point(
            data = grp_sum,
            ggplot2::aes(x = group, y = value),
            colour = col, fill = col,
            shape = 18, size = 5
          )
      } else {
        p <- p +
          ggplot2::geom_violin(
            data = sub,
            fill = col, colour = col, alpha = 0.30,
            trim = TRUE, scale = "width"
          ) +
          .jitter_layer(sub, col)
      }
    }
    p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 1, scales = "free_y")

  } else {
    col <- .pa_metric_cols[unique(d$metric_label)]
    if (.is_const(d)) {
      grp_sum <- aggregate(value ~ group, data = d, FUN = mean)
      p <- p +
        ggplot2::geom_point(
          data = grp_sum,
          ggplot2::aes(x = group, y = value),
          colour = col, fill = col, shape = 18, size = 5
        )
    } else {
      p <- p +
        ggplot2::geom_violin(
          fill = col, colour = col, alpha = 0.30,
          trim = TRUE, scale = "width"
        ) +
        .jitter_layer(d, col)
    }
  }

  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = "Group", y = "Accuracy") +
    theme

  # Add present colour scale when colouring by observed status
  if (has_present_col) {
    p <- p + ggplot2::scale_colour_manual(
      name   = NULL,
      values = .pa_present_cols
    )
  }

  if (return_data) return(list(plot = p, data = d))
  p
}

# =============================================================================
# HEATMAP  (single model, by_variety)
# =============================================================================

.pa_heatmap <- function(res, metric, theme, return_data) {

  d <- .pa_long_bv(res, metric)
  d <- .pa_order_groups(d)
  d <- .pa_order_varieties(d)

  use_facet <- length(unique(d$metric_label)) > 1L

  p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = variety, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(
      option = "D",
      name   = "Accuracy",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    ggplot2::labs(x = "Group", y = "Variety") +
    theme +
    ggplot2::theme(
      axis.text.y      = ggplot2::element_text(size = 7),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Overlay a small grey solid circle on tiles where the variety was not
  # observed at that site (unobserved).
  if ("present" %in% names(d)) {
    d_unobs <- d[!d$present, , drop = FALSE]
    if (nrow(d_unobs) > 0L)
      p <- p + ggplot2::geom_point(
        data        = d_unobs,
        ggplot2::aes(x = group, y = variety),
        colour      = "grey50",
        shape       = 16,
        size        = 1.5,
        inherit.aes = FALSE
      )
  }

  if (use_facet)
    p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2)

  if (return_data) return(list(plot = p, data = d))
  p
}

# =============================================================================
# DUMBBELL  (two models, group-level)
# =============================================================================

.pa_dumbbell <- function(res, res2, metric, label1, label2, theme, return_data) {

  res  <- .pa_bv_to_grp(res)
  res2 <- .pa_bv_to_grp(res2)

  d1         <- .pa_long_grp(res,  metric); d1$model_label <- label1
  d2         <- .pa_long_grp(res2, metric); d2$model_label <- label2
  dall       <- rbind(d1, d2)

  # Segment data: one row per group x metric
  seg <- merge(
    d1[, c("group", "value", "metric_label")],
    d2[, c("group", "value", "metric_label")],
    by = c("group", "metric_label"),
    suffixes = c("_1", "_2")
  )
  names(seg)[names(seg) == "value_1"] <- "x1"
  names(seg)[names(seg) == "value_2"] <- "x2"
  seg$gain_dir <- ifelse(seg$x2 >= seg$x1, "Improved", "Declined")

  # Order groups by model1 accuracy value (ascending -> best at top)
  ref_metric     <- unique(dall$metric_label)[1L]
  ref_d1         <- d1[d1$metric_label == ref_metric, ]
  grp_ord        <- names(sort(tapply(ref_d1$value, ref_d1$group, mean, na.rm = TRUE)))
  dall$group     <- factor(dall$group, levels = grp_ord)
  seg$group      <- factor(seg$group,  levels = grp_ord)
  dall$model_label <- factor(dall$model_label, levels = c(label1, label2))

  use_facet <- length(unique(dall$metric_label)) > 1L

  model_fills <- c("#2166AC", "#D6604D")
  names(model_fills) <- c(label1, label2)

  p <- ggplot2::ggplot() +
    # Connecting segment -- colour encodes direction of change
    ggplot2::geom_segment(
      data = seg,
      ggplot2::aes(x     = x1, xend  = x2,
                   y     = group, yend  = group,
                   colour = gain_dir),
      linewidth = 1.4, alpha = 0.85
    ) +
    # Model points -- shape + fill distinguish models
    ggplot2::geom_point(
      data = dall,
      ggplot2::aes(x    = value,
                   y    = group,
                   shape = model_label,
                   fill  = model_label),
      colour = "white", stroke = 0.5, size = 4.0
    ) +
    ggplot2::scale_colour_manual(
      name   = NULL,
      values = .pa_gain_cols
    ) +
    ggplot2::scale_shape_manual(
      name   = "Model",
      values = c(21L, 22L),
      labels = c(label1, label2)
    ) +
    ggplot2::scale_fill_manual(
      name   = "Model",
      values = model_fills,
      labels = c(label1, label2)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::guides(
      shape  = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          fill   = unname(model_fills),
          colour = "white",
          size   = 3.5
        )
      ),
      fill   = "none",
      colour = ggplot2::guide_legend(
        order = 2,
        override.aes = list(linewidth = 1.4, alpha = 0.85)
      )
    ) +
    ggplot2::labs(x = "Accuracy", y = NULL) +
    theme

  if (use_facet)
    p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2, scales = "free_x")

  if (return_data)
    return(list(plot = p,
                data = list(model1 = d1, model2 = d2, segments = seg)))
  p
}

# =============================================================================
# SCATTER  (two models, by_variety)
# =============================================================================

.pa_scatter <- function(res, res2, metric, label1, label2, theme, return_data) {

  d1 <- .pa_long_bv(res,  metric)
  d2 <- .pa_long_bv(res2, metric)
  names(d1)[names(d1) == "value"] <- ".value_x"
  names(d2)[names(d2) == "value"] <- ".value_y"

  d <- merge(
    d1[, c("group", "variety", "metric_label", ".value_x")],
    d2[, c("group", "variety", "metric_label", ".value_y")],
    by = c("group", "variety", "metric_label")
  )

  use_facet <- length(unique(d$metric_label)) > 1L
  n_grp     <- length(unique(d$group))

  p <- ggplot2::ggplot(d,
       ggplot2::aes(x = .value_x, y = .value_y, colour = group)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", colour = "grey45", linewidth = 0.7
    ) +
    ggplot2::geom_point(alpha = 0.70, size = 2.0) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      x      = paste0(label1, "  accuracy"),
      y      = paste0(label2, "  accuracy"),
      colour = "Group"
    ) +
    theme

  if (use_facet)
    p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2)

  if (return_data) return(list(plot = p, data = d))
  p
}

# =============================================================================
# DIFF  (two models, group-level)
# =============================================================================

.pa_diff <- function(res, res2, metric, label1, label2, theme, return_data) {

  # ---- Single-model present-gap mode ----------------------------------------
  is_bv        <- "variety" %in% names(res)
  has_pres_spl <- is_bv && "present" %in% names(res) &&
                  any(res$present) && any(!res$present)

  if (is.null(res2) && has_pres_spl) {
    d <- .pa_present_grp(res, metric)

    d_obs   <- d[d$present, ,  drop = FALSE]
    d_unobs <- d[!d$present, , drop = FALSE]

    dd <- merge(
      d_obs  [, c("group", "value", "metric_label")],
      d_unobs[, c("group", "value", "metric_label")],
      by = c("group", "metric_label"), suffixes = c("_obs", "_unobs")
    )
    dd$diff_val <- dd$value_unobs - dd$value_obs
    dd$gain_dir <- ifelse(dd$diff_val >= 0, "Improved", "Declined")

    ref_metric <- unique(dd$metric_label)[1L]
    ref_d      <- dd[dd$metric_label == ref_metric, ]
    grp_ord    <- names(sort(tapply(ref_d$diff_val, ref_d$group, mean,
                                    na.rm = TRUE)))
    dd$group   <- factor(dd$group, levels = grp_ord)

    use_facet <- length(unique(dd$metric_label)) > 1L
    abs_max   <- max(abs(dd$diff_val), na.rm = TRUE)
    x_lim     <- c(-abs_max, abs_max) * 1.08

    p <- ggplot2::ggplot(dd, ggplot2::aes(x = diff_val, y = group,
                                           colour = gain_dir)) +
      ggplot2::geom_vline(
        xintercept = 0, linetype = "dashed",
        colour = "grey40", linewidth = 0.6
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = diff_val, yend = group),
        linewidth = 0.9
      ) +
      ggplot2::geom_point(size = 3.5) +
      ggplot2::scale_colour_manual(values = .pa_gain_cols, guide = "none") +
      ggplot2::scale_x_continuous(
        limits = x_lim,
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::labs(
        x = "Accuracy gap  (Unobserved \u2212 Observed)",
        y = NULL
      ) +
      theme

    if (use_facet)
      p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2, scales = "free_x")

    if (return_data) return(list(plot = p, data = dd))
    return(p)
  }

  # ---- Two-model comparison mode --------------------------------------------
  res  <- .pa_bv_to_grp(res)
  res2 <- .pa_bv_to_grp(res2)
  d1 <- .pa_long_grp(res,  metric)
  d2 <- .pa_long_grp(res2, metric)

  d <- merge(
    d1[, c("group", "value", "n_vars", "metric_label")],
    d2[, c("group", "value", "metric_label")],
    by       = c("group", "metric_label"),
    suffixes = c("_1", "_2")
  )
  d$diff_val <- d$value_2 - d$value_1
  d$gain_dir <- ifelse(d$diff_val >= 0, "Improved", "Declined")

  # Order by diff_val within the first metric (ascending -> biggest gain at top)
  ref_metric <- unique(d$metric_label)[1L]
  ref_d      <- d[d$metric_label == ref_metric, ]
  grp_ord    <- names(sort(tapply(ref_d$diff_val, ref_d$group, mean, na.rm = TRUE)))
  d$group    <- factor(d$group, levels = grp_ord)

  use_facet <- length(unique(d$metric_label)) > 1L

  # Symmetric x limits centred on 0
  abs_max <- max(abs(d$diff_val), na.rm = TRUE)
  x_lim   <- c(-abs_max, abs_max) * 1.08

  p <- ggplot2::ggplot(d, ggplot2::aes(x = diff_val, y = group,
                                        colour = gain_dir)) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      colour = "grey40", linewidth = 0.6
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = diff_val, yend = group),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_manual(values = .pa_gain_cols, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = x_lim,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      x = sprintf("Accuracy gain  (%s \u2212 %s)", label2, label1),
      y = NULL
    ) +
    theme

  if (use_facet)
    p <- p + ggplot2::facet_wrap(~ metric_label, ncol = 2, scales = "free_x")

  if (return_data) return(list(plot = p, data = d))
  p
}
