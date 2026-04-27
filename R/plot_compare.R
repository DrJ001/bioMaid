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
## hover_data is an optional list of per-panel metadata for hover interaction.
#' @noRd
.pc_make_interactive <- function(p, hover_data = NULL) {
  structure(list(plot = p, hover_data = hover_data), class = "pc_interactive")
}

## Support ggplot2's + operator on pc_interactive objects.
##
## ggplot2 >= 4.0.0 migrated to S7 classes and registers an S7 method for
## the base + operator that returns NULL for non-ggplot left operands.
## Neither +.pc_interactive nor Ops.pc_interactive can intercept that S7
## dispatch when the right operand is an S7 gg object.
##
## The fix: use ggplot2::ggplot_add() directly instead of the + operator
## on the internal ggplot.  ggplot_add(object, plot, name) is the exported
## ggplot2 S7 generic that + calls internally, so the effect is identical
## but the call goes through ggplot2's own API without hitting the dispatch
## conflict.
#' @export
Ops.pc_interactive <- function(e1, e2) {
  if (.Generic == "+" && inherits(e1, "pc_interactive")) {
    objectname <- deparse(substitute(e2))
    e1$plot <- ggplot2::ggplot_add(e2, e1$plot, objectname)
    return(e1)
  }
  stop("Operator '", .Generic, "' is not defined for pc_interactive objects.")
}

## Auto-convert to plotly on print.
## For dotplots, injects hover-band JavaScript if hover_data is present.
#' @export
print.pc_interactive <- function(x, ...) {
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required to display an interactive plot. ",
         "Install with: install.packages('plotly')")
  p_plotly <- plotly::ggplotly(x$plot, tooltip = "text")
  if (!is.null(x$hover_data))
    p_plotly <- .pc_add_hover_js(p_plotly, x$hover_data)
  print(p_plotly)
  invisible(x)
}

## knitr knit_print method for pc_interactive.
##
## When a pc_interactive object is the last expression in a knitr/Quarto
## chunk, knitr dispatches to knit_print() rather than print().  Without
## this method knitr would call print.pc_interactive() which prints the
## plotly widget via a nested print() call — a side-effect that knitr
## cannot capture — and the chunk output shows NULL instead of the widget.
## This method converts to plotly and hands the widget back to knitr via
## knitr::knit_print() so it is properly embedded in the HTML output.
##
## @importFrom knitr knit_print is required so roxygen2 registers this as
## S3method(knit_print, pc_interactive) in NAMESPACE rather than a plain
## export — without the S3method entry, UseMethod("knit_print") in knitr
## will never dispatch here.
#' @importFrom knitr knit_print
#' @export
knit_print.pc_interactive <- function(x, ...) {
  if (!requireNamespace("plotly",  quietly = TRUE))
    stop("Package 'plotly' is required to display an interactive plot. ",
         "Install with: install.packages('plotly')")
  if (!requireNamespace("knitr",   quietly = TRUE))
    stop("Package 'knitr' is required for knit_print support.")
  p_plotly <- plotly::ggplotly(x$plot, tooltip = "text")
  if (!is.null(x$hover_data))
    p_plotly <- .pc_add_hover_js(p_plotly, x$hover_data)
  knitr::knit_print(p_plotly, ...)
}

## Convert a pc_interactive object to a plotly htmlwidget.
##
## as_plotly() is the escape hatch for knitr/Quarto vignettes where S3
## dispatch through knit_print may be unreliable (e.g. when the installed
## package version lags the source).  Calling as_plotly() as the LAST
## expression of a chunk returns a native plotly htmlwidget, which knitr
## captures via htmlwidgets' built-in knit_print.htmlwidget — no biomAid-
## specific registration required.
##
## Usage in a vignette chunk:
##   p_int <- plot_compare(res, interactive = TRUE) + ggtitle("...")
##   as_plotly(p_int)   # last expression: knitr embeds the widget directly

#' Convert a `pc_interactive` Object to a Plotly Widget
#'
#' @description
#' Converts a [plot_compare()] `pc_interactive` result to a
#' [plotly::plotly] htmlwidget, applying the hover-band JavaScript for
#' dotplots.  Useful in knitr/Quarto documents where the object should be
#' the direct return value of the chunk so that knitr embeds it correctly.
#'
#' @param x A `pc_interactive` object returned by
#'   `plot_compare(..., interactive = TRUE)`.
#' @param height Numeric or `NULL`. Widget height in pixels passed to
#'   [plotly::ggplotly()].  `NULL` (default) lets plotly choose its own
#'   default height (~400 px), which is often too short when the plot has
#'   multiple facet panels each showing many varieties.  A sensible rule of
#'   thumb: 180 px per panel row (e.g. `height = 700` for 4 panels in a
#'   2 × 2 grid).
#' @param width Numeric or `NULL`. Widget width in pixels.  `NULL` (default)
#'   uses 100 % of the container width.
#' @param ... Unused; included for future extensibility.
#'
#' @return A [plotly::plotly] htmlwidget.
#'
#' @examples
#' \dontrun{
#' # Single-group dotplot — default height is fine
#' as_plotly(plot_compare(res, type = "dotplot", interactive = TRUE))
#'
#' # Four-panel dotplot — set a taller height so panels don't overlap
#' as_plotly(plot_compare(res, type = "dotplot", interactive = TRUE),
#'           height = 700)
#' }
#'
#' @export
as_plotly <- function(x, height = NULL, width = NULL, ...) {
  if (!inherits(x, "pc_interactive"))
    stop("'x' must be a pc_interactive object returned by plot_compare().")
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  p_plotly <- plotly::ggplotly(x$plot, tooltip = "text",
                               height = height, width = width)
  p_plotly <- .pc_fix_subplot_spacing(p_plotly)
  if (!is.null(x$hover_data))
    p_plotly <- .pc_add_hover_js(p_plotly, x$hover_data)
  p_plotly
}

## Fix subplot y-axis domain spacing for multi-panel ggplotly output.
##
## ggplotly converts facet_wrap panels into plotly subplots, but the computed
## domain gaps are often too small — panels crowd each other regardless of
## the total widget height.  This function redistributes yaxis domains so
## that each row occupies an equal share of [0, 1] with a fixed fractional
## gap between rows, then shifts any strip-label annotations (yref="paper")
## that sat at the top of each old domain to the top of the new domain.
#' @noRd
.pc_fix_subplot_spacing <- function(p_plotly, row_gap = 0.12) {

  ly <- p_plotly$x$layout

  # ---- Identify multi-row yaxis entries ----------------------------------
  yax_names <- sort(names(ly)[grepl("^yaxis[0-9]*$", names(ly))])
  if (length(yax_names) <= 1L) return(p_plotly)   # single panel — nothing to do

  old_domains <- lapply(yax_names, function(n) ly[[n]]$domain)
  old_mids    <- vapply(old_domains, function(d)
                          if (is.null(d)) NA_real_ else mean(d), numeric(1L))

  # Cluster into rows: round mid to 1 d.p. then rank ascending (1=bottom row)
  row_id <- as.integer(factor(round(old_mids, 1L)))
  n_rows <- max(row_id, na.rm = TRUE)
  if (n_rows <= 1L) return(p_plotly)              # all panels in one row

  # ---- Compute new evenly-spaced row domains -----------------------------
  row_h  <- (1 - row_gap * (n_rows - 1L)) / n_rows
  row_lo <- function(r) (r - 1L) * (row_h + row_gap)
  row_hi <- function(r) row_lo(r) + row_h

  # Map from old domain top → new domain top (used to relocate annotations)
  old_top_to_new <- stats::setNames(
    vapply(row_id, row_hi, numeric(1L)),
    vapply(old_domains, function(d) if (is.null(d)) NA_real_ else d[2L],
           numeric(1L))
  )

  # Apply new domains to yaxis entries
  for (i in seq_along(yax_names)) {
    r  <- row_id[i]
    p_plotly$x$layout[[yax_names[i]]]$domain <- c(row_lo(r), row_hi(r))
  }

  # ---- Relocate strip-label annotations ----------------------------------
  # ggplotly places strip text as annotations with yref="paper" whose y
  # coordinate sits just above the top of the panel's old yaxis domain.
  # Shift them to the top of the new domain (+/- the same small offset).
  anns <- p_plotly$x$layout$annotations
  if (!is.null(anns)) {
    tol <- 0.08   # tolerance: annotation is "above" a domain top if within tol
    p_plotly$x$layout$annotations <- lapply(anns, function(a) {
      if (!identical(a$yref, "paper")) return(a)
      ay <- a$y
      if (is.null(ay) || is.na(ay)) return(a)
      # Find the old domain top this annotation belongs to
      match_top <- names(old_top_to_new)[
        abs(as.numeric(names(old_top_to_new)) - ay) < tol
      ]
      if (length(match_top) == 1L) {
        old_top <- as.numeric(match_top)
        new_top <- old_top_to_new[[match_top]]
        a$y <- new_top + (ay - old_top)   # preserve offset above panel top
      }
      a
    })
  }

  p_plotly
}

## Inject hover-band JavaScript into a plotly dotplot.
##
## The geom_rect band is deliberately excluded from the ggplot when building
## for interactive mode (the builder receives interactive = TRUE).  This means
## the ggplotly output contains NO fill="toself" band traces to remove, so no
## ggplot theming elements are accidentally stripped.  The band is added here
## as a proper layout.shape which Plotly.relayout() CAN update at runtime.
##
## hover_data: named list, keys = plotly xaxis names ("x", "x2", ...)
##   Each element has: crit, orig_x0, orig_x1, n_vars
#' @noRd
.pc_add_hover_js <- function(p_plotly, hover_data) {

  # ---- APPEND band shapes to existing layout.shapes ----------------------
  # ggplotly already puts shapes in layout.shapes (panel borders, strip
  # boxes, etc.).  We must APPEND the band shapes — not replace — or all
  # those theming shapes are lost.
  existing <- p_plotly$x$layout$shapes
  if (is.null(existing)) existing <- list()
  n_existing <- length(existing)          # 0-based offset for our new shapes

  xax_names <- names(hover_data)
  band_shapes <- lapply(seq_along(xax_names), function(i) {
    xax   <- xax_names[i]
    yax   <- if (xax == "x") "y" else paste0("y", sub("^x", "", xax))
    panel <- hover_data[[xax]]
    list(
      type      = "rect",
      xref      = xax,
      yref      = yax,
      x0        = panel$orig_x0,
      x1        = panel$orig_x1,
      y0        = 0.5,
      y1        = panel$n_vars + 0.5,
      fillcolor = "rgba(255,239,192,0.6)",
      line      = list(width = 0L),
      layer     = "below"
    )
  })
  p_plotly$x$layout$shapes <- c(existing, band_shapes)

  # Build an explicit xaxis -> band shape index (0-based) map so the
  # JavaScript never has to scan for the right shape by type — it knows
  # exactly which indices are the bands we just appended.
  axis_shape_map <- stats::setNames(
    as.list(seq_len(length(xax_names)) + n_existing - 1L),
    xax_names
  )

  # ---- Inject hover JavaScript -------------------------------------------
  panels_json    <- jsonlite::toJSON(hover_data,      auto_unbox = TRUE)
  axis_shape_json <- jsonlite::toJSON(axis_shape_map, auto_unbox = TRUE)

  js <- paste0("
function(el, x) {

  // Per-panel data: criterion + original band bounds (embedded from R)
  var panelData = ", panels_json, ";

  // Explicit band shape indices (0-based) per panel axis — embedded from R
  // so the JS never accidentally selects a ggplotly theming shape.
  var axisShape = ", axis_shape_json, ";

  // Map xaxis -> [traceIdx, ...] for scatter+marker traces.
  // Exclude reference diamonds and line-only traces.
  var axisTraces = {};
  var origColours = {};
  x.data.forEach(function(trace, idx) {
    if (trace.type !== 'scatter') return;
    var mode = trace.mode || '';
    if (mode.indexOf('markers') < 0) return;        // skip line-only traces
    var sym = trace.marker ? trace.marker.symbol : null;
    if (sym === 'diamond' || sym === 18) return;    // skip reference diamond
    var xax = trace.xaxis || 'x';
    if (!axisTraces[xax]) axisTraces[xax] = [];
    axisTraces[xax].push(idx);
    origColours[idx] = (trace.marker && trace.marker.color !== undefined)
                         ? trace.marker.color : '#4E79A7';
  });

  // ---- Hover: move band + recolour points --------------------------------
  el.on('plotly_hover', function(eventData) {
    if (!eventData.points || !eventData.points.length) return;
    var pt      = eventData.points[0];
    var traceEl = el._fullData[pt.curveNumber];
    var xax     = (traceEl && traceEl.xaxis) ? traceEl.xaxis : 'x';
    var panel   = panelData[xax];
    if (!panel) return;

    var hX   = pt.x;
    var crit = panel.crit;

    // Move the band shape
    var si = axisShape[xax];
    if (si !== undefined) {
      var su = {};
      su['shapes[' + si + '].x0'] = hX - crit;
      su['shapes[' + si + '].x1'] = hX + crit;
      Plotly.relayout(el, su);
    }

    // Recolour all marker points in this panel
    (axisTraces[xax] || []).forEach(function(ti) {
      var td = el._fullData[ti];
      if (!td || !td.x) return;
      var cols = td.x.map(function(pv) {
        if (pv > hX + crit) return '#59A14F';   // sig. better (green)
        if (pv < hX - crit) return '#E15759';   // sig. worse  (red)
        return '#4E79A7';                         // not sig.    (blue)
      });
      Plotly.restyle(el, {'marker.color': [cols]}, [ti]);
    });
  });

  // ---- Unhover: restore original state -----------------------------------
  el.on('plotly_unhover', function() {
    // Restore each panel's band
    Object.keys(axisShape).forEach(function(xax) {
      var panel = panelData[xax];
      if (!panel) return;
      var si = axisShape[xax];
      var su = {};
      su['shapes[' + si + '].x0'] = panel.orig_x0;
      su['shapes[' + si + '].x1'] = panel.orig_x1;
      Plotly.relayout(el, su);
    });
    // Restore each trace's original colours
    Object.keys(origColours).forEach(function(ti) {
      var col = origColours[ti];
      Plotly.restyle(el, {'marker.color': col}, [parseInt(ti)]);
    });
  });
}
")

  htmlwidgets::onRender(p_plotly, js)
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
## Two-pass strategy:
##
## Pass 1 — criterion-based grouping.
##   Group rows by their criterion value.  Columns that are constant within
##   every criterion group are by-columns; the one that varies is the
##   comparison variable.  This gives the correct answer in the common case
##   where each by-group has a different criterion (e.g. 4 sites, 4 HSD
##   values → Treatment varies within each group, Site is constant).
##
## Pass 2 — cardinality fallback.
##   Used only when Pass 1 returns more than one comparison candidate
##   (i.e. the criterion is identical across all groups).  Among the
##   candidates that varied in Pass 1, the one with the most unique values
##   is the comparison variable.
##
##   Example where Pass 1 would fail: 2 sites with the same HSD = 180 → one
##   big criterion group → both Site and Variety appear to vary → cardinality
##   correctly picks Variety (20 unique) over Site (2 unique).
#' @noRd
.pc_detect_cols <- function(res, crit_col) {

  skip       <- c("predicted.value", "std.error", "avsed", crit_col)
  candidates <- setdiff(names(res), skip)

  if (length(candidates) == 0L)
    stop("No factor columns found in 'res' to use as the comparison variable.")

  # ---- Pass 1: criterion-based grouping ----------------------------------
  grp_id    <- as.character(res[[crit_col]])
  by_crit   <- character(0L)
  comp_crit <- character(0L)

  for (col in candidates) {
    vals_per_grp <- tapply(as.character(res[[col]]), grp_id,
                           function(x) length(unique(x)))
    if (all(vals_per_grp == 1L))
      by_crit   <- c(by_crit,   col)
    else
      comp_crit <- c(comp_crit, col)
  }

  if (length(comp_crit) == 1L)
    return(list(by_cols = by_crit, comp_col = comp_crit))

  # ---- Pass 2: cardinality tiebreaker ------------------------------------
  # Restrict to the candidates that varied in Pass 1 (or all if Pass 1 found
  # nothing, which should not happen in practice).
  pool <- if (length(comp_crit) > 0L) comp_crit else candidates

  n_unique <- vapply(pool,
                     function(col) length(unique(as.character(res[[col]]))),
                     integer(1L))

  max_u    <- max(n_unique)
  comp_col <- pool[n_unique == max_u]

  if (length(comp_col) > 1L) {
    warning("Could not unambiguously identify the comparison variable. ",
            "Candidates: ", paste(comp_col, collapse = ", "),
            ". Using '", comp_col[1L], "'.", call. = FALSE)
    comp_col <- comp_col[1L]
  }

  by_cols <- setdiff(candidates, comp_col)

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

  res <- res[order(res$Group, res$y_pos), ]

  # var_label: group-namespaced ordered factor for static discrete y axis.
  # Format "Group\x01Variety" (SOH separator) lets scale_y_discrete strip
  # the prefix in labels.  \x01 (ASCII SOH) is safe in R string literals
  # and cannot appear in group or variety names produced by biomAid.
  # Levels are ordered by y_pos ascending within each group (y_pos 1 = lowest
  # pred = bottom of panel, y_pos n = highest pred = top of panel).
  # This is the ONLY correct way to get per-panel variety labels with
  # facet_wrap(scales="free"): unlist(grp_labels) prepends the list element
  # name, so bare position lookups return NA for all panels.
  res$var_label <- factor(
    paste0(res$Group, "\x01", res$Variety),
    levels = unique(paste0(res$Group, "\x01", res$Variety))
  )

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
.pc_plot_dotplot <- function(df, crit_col, reference, n_groups,
                              interactive = FALSE, theme, ...) {

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

  # ---- y-axis strategy: static vs interactive --------------------------------
  #
  # STATIC (interactive = FALSE):
  #   Use var_label (a "Group\x00Variety" ordered factor) as the y aesthetic.
  #   scale_y_discrete strips the group prefix in labels.  Each facet panel
  #   automatically shows only its own group's factor levels, in y_pos order.
  #   geom_rect uses ymin = -Inf / ymax = Inf to span the full panel height,
  #   which works on a discrete scale and avoids passing finite band bounds.
  #
  # INTERACTIVE (interactive = TRUE):
  #   Keep y = y_pos (numeric) so band layout.shapes can be positioned and
  #   moved by Plotly.relayout() in data coordinates.  The label function
  #   (grp_labels[[1L]] lookup) is correct for single-group plots; for
  #   multi-group interactive the hover tooltip still shows variety names.

  if (!interactive) {

    p <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = var_label))

    # Criterion band: spans full panel height on the discrete scale
    p <- p + ggplot2::geom_rect(
      data        = band_df,
      ggplot2::aes(xmin = band_lo, xmax = band_hi,
                   ymin = -Inf,    ymax = Inf),
      fill        = "#FFEFC0",
      alpha       = 0.6,
      inherit.aes = FALSE
    )

    p <- p +
      ggplot2::geom_vline(
        data      = band_df,
        ggplot2::aes(xintercept = ref_pred),
        linetype  = "dashed",
        linewidth = 0.45,
        colour    = "grey40"
      ) +
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
      ggplot2::scale_y_discrete(
        # Strip the "Group\x01" namespace prefix — show only the variety name
        labels  = function(x) sub("^[^\x01]+\x01", "", x),
        expand  = ggplot2::expansion(add = 0.6)
      ) +
      ggplot2::labs(x = "Predicted value", y = NULL, caption = caption) +
      theme +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey92",
                                                 colour = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold"),
        legend.position  = "bottom",
        panel.spacing    = ggplot2::unit(0.5, "lines"),
        axis.text.y      = ggplot2::element_text(size = 7)
      )

    # Reference variety: diamond overlay using var_label for y
    if (using_ref) {
      ref_rows <- df[df$Variety == reference, ]
      if (nrow(ref_rows) > 0L) {
        ref_rows$tooltip_text <- paste0(
          "<b>", ref_rows$Variety, "</b> (reference)",
          "<br>Predicted: ", round(ref_rows$pred, 2L),
          "<br>Group: ",     ref_rows$Group
        )
        p <- p +
          ggplot2::geom_point(
            data        = ref_rows,
            ggplot2::aes(x = pred, y = var_label, text = tooltip_text),
            shape       = 18L,
            size        = 4.5,
            colour      = "grey20",
            inherit.aes = FALSE
          )
      }
    }

  } else {

    # --- Interactive branch: numeric y_pos --------------------------------
    grp_labels <- lapply(unique(df$Group), function(g) {
      sub <- df[df$Group == g, ]
      stats::setNames(sub$Variety, sub$y_pos)
    })
    names(grp_labels) <- unique(df$Group)
    all_y    <- sort(unique(df$y_pos))
    label_fn <- function(breaks) {
      lbl <- grp_labels[[1L]][as.character(breaks)]
      ifelse(is.na(lbl), "", lbl)
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = pred, y = y_pos))

    p <- p +
      ggplot2::geom_vline(
        data      = band_df,
        ggplot2::aes(xintercept = ref_pred),
        linetype  = "dashed",
        linewidth = 0.45,
        colour    = "grey40"
      ) +
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
      ggplot2::labs(x = "Predicted value", y = NULL, caption = caption) +
      theme +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey92",
                                                 colour = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold"),
        legend.position  = "bottom",
        panel.spacing    = ggplot2::unit(0.5, "lines"),
        axis.text.y      = ggplot2::element_text(size = 7)
      )

    # Reference variety: diamond overlay using y_pos
    if (using_ref) {
      ref_rows <- df[df$Variety == reference, ]
      if (nrow(ref_rows) > 0L) {
        ref_rows$tooltip_text <- paste0(
          "<b>", ref_rows$Variety, "</b> (reference)",
          "<br>Predicted: ", round(ref_rows$pred, 2L),
          "<br>Group: ",     ref_rows$Group
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

  }  # end static / interactive branch

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

  # Use var_label (discrete, group-namespaced) for the y axis so that
  # each facet panel shows its own correctly ordered variety labels.
  # The scale_y_discrete label function strips the "Group\x00" prefix.
  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = pred, y = var_label)) +
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
    ggplot2::scale_y_discrete(
      labels = function(x) sub("^[^\x01]+\x01", "", x),
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
    dotplot = .pc_plot_dotplot(df, crit_col, reference, n_groups,
                               interactive, theme, ...),
    letters = .pc_plot_letters(df, crit_col, n_groups, theme, ...),
    heatmap = .pc_plot_heatmap(df, crit_col, n_groups, theme, ...)
  ))

  # ---- Optionally wrap for lazy plotly conversion -------------------------
  # Return a pc_interactive wrapper so the user can still use + to add
  # ggplot2 layers before the plotly conversion happens at print() time.
  if (interactive) {
    # For dotplots, build per-panel hover metadata so the JavaScript can
    # move the criterion band and recolour points on hover.
    hover_data <- if (type == "dotplot") {
      grps <- unique(df$Group)
      pdata <- lapply(seq_along(grps), function(gi) {
        g   <- grps[gi]
        sub <- df[df$Group == g, ]
        list(
          crit    = as.numeric(unique(sub$crit)[1L]),
          orig_x0 = as.numeric(unique(sub$band_lo)[1L]),
          orig_x1 = as.numeric(unique(sub$band_hi)[1L]),
          n_vars  = as.integer(nrow(sub))   # needed for shape y1 bound
        )
      })
      # Plotly names axes "x", "x2", "x3", ...
      xax_names <- ifelse(seq_along(grps) == 1L, "x",
                          paste0("x", seq_along(grps)))
      stats::setNames(pdata, xax_names)
    } else {
      NULL
    }
    return(.pc_make_interactive(p, hover_data))
  }

  p
}
