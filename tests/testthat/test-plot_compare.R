# tests/testthat/test-plot_compare.R
# Tests for plot_compare() — no ASReml needed; uses synthetic compare() output.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Build a synthetic compare() result for a single group (no by).
make_compare_single <- function(n = 10L, type = "HSD", seed = 1L) {
  set.seed(seed)
  preds <- sort(rnorm(n, 50, 8), decreasing = FALSE)
  df <- data.frame(
    Variety         = factor(paste0("G", sprintf("%02d", seq_len(n)))),
    predicted.value = preds,
    std.error       = runif(n, 0.5, 1.5),
    stringsAsFactors = FALSE
  )
  df[[type]]  <- rep(10.0, n)   # constant criterion
  df$avsed    <- rep(5.0,  n)
  df
}

# Build a synthetic compare() result with a 'by' grouping.
make_compare_grouped <- function(groups = c("S1", "S2"),
                                 n_var  = 8L,
                                 type   = "HSD",
                                 seed   = 2L) {
  set.seed(seed)
  rows <- lapply(groups, function(g) {
    preds <- sort(rnorm(n_var, 50, 8), decreasing = FALSE)
    df <- data.frame(
      Site            = g,
      Variety         = paste0("V", sprintf("%02d", seq_len(n_var))),
      predicted.value = preds,
      std.error       = runif(n_var, 0.5, 1.5),
      stringsAsFactors = FALSE
    )
    df[[type]] <- rep(runif(1L, 8, 14), n_var)
    df$avsed   <- rep(runif(1L, 4,  7),  n_var)
    df
  })
  out       <- do.call(rbind, rows)
  out$Site  <- factor(out$Site)
  out
}

# Build a two-factor grouped result (Site x Treatment as by).
make_compare_2factor <- function(type = "HSD", seed = 3L) {
  set.seed(seed)
  combos <- expand.grid(Site = c("S1", "S2"),
                        Treatment = c("T0", "T1"),
                        stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(combos)), function(i) {
    g     <- combos[i, ]
    n_var <- 6L
    preds <- sort(rnorm(n_var, 50, 8))
    df    <- data.frame(
      Site            = g$Site,
      Treatment       = g$Treatment,
      Variety         = paste0("V", sprintf("%02d", seq_len(n_var))),
      predicted.value = preds,
      std.error       = runif(n_var, 0.5, 1.5),
      stringsAsFactors = FALSE
    )
    df[[type]] <- rep(runif(1L, 8, 14), n_var)
    df$avsed   <- rep(runif(1L, 4,  7),  n_var)
    df
  })
  out           <- do.call(rbind, rows)
  out$Site      <- factor(out$Site)
  out$Treatment <- factor(out$Treatment)
  out
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("non-data-frame input stops", {
  expect_error(plot_compare(list()), "data frame")
})

test_that("missing required columns stops", {
  df <- data.frame(Variety = "G1", HSD = 1)
  expect_error(plot_compare(df), "missing expected column")
})

test_that("no criterion column stops", {
  res <- make_compare_single()
  res$HSD <- NULL
  expect_error(plot_compare(res), "No criterion column")
})

test_that("invalid type stops", {
  res <- make_compare_single()
  expect_error(plot_compare(res, type = "bar"), "should be one of")
})

test_that("non-theme 'theme' argument stops", {
  res <- make_compare_single()
  expect_error(plot_compare(res, theme = "bw"), "'theme' must be a ggplot2")
})

test_that("non-logical return_data stops", {
  res <- make_compare_single()
  expect_error(plot_compare(res, return_data = "yes"), "single logical")
})

# ---------------------------------------------------------------------------
# 2. Return types — single group
# ---------------------------------------------------------------------------

test_that("dotplot returns a ggplot object (single group)", {
  res <- make_compare_single()
  p   <- plot_compare(res, type = "dotplot")
  expect_s3_class(p, "ggplot")
})

test_that("letters returns a ggplot object (single group)", {
  res <- make_compare_single()
  p   <- plot_compare(res, type = "letters")
  expect_s3_class(p, "ggplot")
})

test_that("heatmap returns a ggplot object (single group)", {
  res <- make_compare_single()
  p   <- plot_compare(res, type = "heatmap")
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# 3. Return types — grouped
# ---------------------------------------------------------------------------

test_that("dotplot returns ggplot for grouped result", {
  res <- make_compare_grouped()
  p   <- plot_compare(res, type = "dotplot")
  expect_s3_class(p, "ggplot")
})

test_that("letters returns ggplot for grouped result", {
  res <- make_compare_grouped()
  p   <- plot_compare(res, type = "letters")
  expect_s3_class(p, "ggplot")
})

test_that("heatmap returns ggplot for grouped result", {
  res <- make_compare_grouped()
  p   <- plot_compare(res, type = "heatmap")
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# 4. return_data = TRUE
# ---------------------------------------------------------------------------

test_that("return_data=TRUE returns a data frame for dotplot", {
  res <- make_compare_single()
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  expect_s3_class(df, "data.frame")
})

test_that("return_data=TRUE returns a data frame for letters", {
  res <- make_compare_single()
  df  <- plot_compare(res, type = "letters", return_data = TRUE)
  expect_s3_class(df, "data.frame")
})

test_that("return_data=TRUE returns a data frame for heatmap", {
  res <- make_compare_single()
  df  <- plot_compare(res, type = "heatmap", return_data = TRUE)
  expect_s3_class(df, "data.frame")
})

test_that("dotplot data has expected columns", {
  res <- make_compare_single()
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  expect_true(all(c("Group", "Variety", "pred", "crit", "rank",
                    "ref_pred", "band_lo", "band_hi", "sig") %in% names(df)))
})

test_that("letters data has letter column", {
  res <- make_compare_single()
  df  <- plot_compare(res, type = "letters", return_data = TRUE)
  expect_true("letter" %in% names(df))
})

test_that("heatmap data has Var_x, Var_y, diff, sig columns", {
  res <- make_compare_single()
  df  <- plot_compare(res, type = "heatmap", return_data = TRUE)
  expect_true(all(c("Var_x", "Var_y", "diff", "sig") %in% names(df)))
})

# ---------------------------------------------------------------------------
# 5. Dotplot logic
# ---------------------------------------------------------------------------

test_that("dotplot: rank 1 is the highest predicted value", {
  res <- make_compare_single(n = 10L, seed = 5L)
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  top <- df[df$rank == 1L, ]
  expect_equal(top$pred, max(df$pred))
})

test_that("dotplot: ref_pred equals max predicted value in group (default)", {
  res <- make_compare_single(n = 10L, seed = 6L)
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  expect_equal(unique(df$ref_pred), max(df$pred))
})

test_that("dotplot: sig='ns' for top-ranked variety (default)", {
  res <- make_compare_single(n = 10L, seed = 7L)
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  expect_equal(df$sig[df$rank == 1L], "ns")
})

test_that("dotplot: varieties below criterion are sig='sig' (default)", {
  # Make criterion small so bottom variety is always significant
  res       <- make_compare_single(n = 10L, seed = 8L)
  res$HSD   <- rep(0.01, 10L)   # tiny criterion
  df        <- plot_compare(res, type = "dotplot", return_data = TRUE)
  # All but rank-1 should be significant
  non_top <- df[df$rank > 1L, ]
  expect_true(all(non_top$sig == "sig"))
})

test_that("dotplot: with large criterion all varieties are sig='ns' (default)", {
  res       <- make_compare_single(n = 10L, seed = 9L)
  res$HSD   <- rep(1000, 10L)
  df        <- plot_compare(res, type = "dotplot", return_data = TRUE)
  expect_true(all(df$sig == "ns"))
})

# ---------------------------------------------------------------------------
# 6. Letters logic
# ---------------------------------------------------------------------------

test_that("letters: all varieties get a non-empty letter string", {
  res <- make_compare_single(n = 8L, seed = 10L)
  df  <- plot_compare(res, type = "letters", return_data = TRUE)
  expect_true(all(nchar(df$letter) > 0L))
})

test_that("letters: single variety gets letter 'a'", {
  lmap <- biomAid:::.pc_cld_group(50, "G1", 10)
  expect_equal(unname(lmap), "a")
})

test_that("letters: when criterion is huge all varieties share letter 'a'", {
  preds <- c(10, 20, 30, 40, 50)
  lmap  <- biomAid:::.pc_cld_group(preds, paste0("G", seq_along(preds)), 1000)
  expect_true(all(lmap == "a"))
})

test_that("letters: when criterion is tiny top and bottom differ", {
  preds <- c(10, 50)
  lmap  <- biomAid:::.pc_cld_group(preds, c("low", "high"), 0.001)
  expect_false(lmap["low"] == lmap["high"])
})

test_that("letters: letter column has same length as input rows", {
  res <- make_compare_single(n = 12L, seed = 11L)
  df  <- plot_compare(res, type = "letters", return_data = TRUE)
  expect_equal(nrow(df), 12L)
})

# ---------------------------------------------------------------------------
# 7. Heatmap logic
# ---------------------------------------------------------------------------

test_that("heatmap data has n^2 rows for single group of n varieties", {
  n   <- 8L
  res <- make_compare_single(n = n, seed = 12L)
  df  <- plot_compare(res, type = "heatmap", return_data = TRUE)
  expect_equal(nrow(df), n^2L)
})

test_that("heatmap: diagonal differences are zero", {
  res  <- make_compare_single(n = 6L, seed = 13L)
  df   <- plot_compare(res, type = "heatmap", return_data = TRUE)
  diag_rows <- as.character(df$Var_x) == as.character(df$Var_y)
  expect_true(all(df$diff[diag_rows] == 0))
})

test_that("heatmap: diagonal entries are never significant", {
  res  <- make_compare_single(n = 6L, seed = 14L)
  df   <- plot_compare(res, type = "heatmap", return_data = TRUE)
  diag_rows <- as.character(df$Var_x) == as.character(df$Var_y)
  expect_true(all(!df$sig[diag_rows]))
})

test_that("heatmap: diff is symmetric", {
  n   <- 5L
  res <- make_compare_single(n = n, seed = 15L)
  df  <- plot_compare(res, type = "heatmap", return_data = TRUE)
  for (v1 in levels(df$Var_x)) {
    for (v2 in levels(df$Var_x)) {
      d_xy <- df$diff[as.character(df$Var_x) == v1 & as.character(df$Var_y) == v2]
      d_yx <- df$diff[as.character(df$Var_x) == v2 & as.character(df$Var_y) == v1]
      if (length(d_xy) && length(d_yx))
        expect_equal(d_xy, d_yx)
    }
  }
})

test_that("heatmap: with tiny criterion all off-diagonal pairs are significant", {
  res     <- make_compare_single(n = 5L, seed = 16L)
  res$HSD <- rep(0.001, 5L)
  df      <- plot_compare(res, type = "heatmap", return_data = TRUE)
  off_diag <- as.character(df$Var_x) != as.character(df$Var_y)
  expect_true(all(df$sig[off_diag]))
})

# ---------------------------------------------------------------------------
# 8. Grouped / faceting
# ---------------------------------------------------------------------------

test_that("dotplot data has correct Group values for 2-group result", {
  res <- make_compare_grouped(groups = c("S1", "S2"))
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  expect_setequal(unique(df$Group), c("S1", "S2"))
})

test_that("heatmap data has correct Group values for 2-group result", {
  res <- make_compare_grouped(groups = c("S1", "S2"))
  df  <- plot_compare(res, type = "heatmap", return_data = TRUE)
  expect_setequal(unique(df$Group), c("S1", "S2"))
})

test_that("2-factor by produces composite group labels in dotplot data", {
  res <- make_compare_2factor()
  df  <- plot_compare(res, type = "dotplot", return_data = TRUE)
  # Expect labels like "S1:T0"
  expect_true(any(grepl(":", df$Group)))
})

test_that("2-factor by: 4 groups present in heatmap data", {
  res <- make_compare_2factor()
  df  <- plot_compare(res, type = "heatmap", return_data = TRUE)
  expect_equal(length(unique(df$Group)), 4L)
})

# ---------------------------------------------------------------------------
# 9. reference argument
# ---------------------------------------------------------------------------

test_that("reference: valid variety name changes ref_pred to that variety's pred", {
  res <- make_compare_single(n = 10L, seed = 30L)
  ref_var  <- as.character(res$Variety[5L])
  ref_pred <- res$predicted.value[5L]
  df  <- plot_compare(res, type = "dotplot", reference = ref_var,
                      return_data = TRUE)
  expect_equal(unique(df$ref_pred), ref_pred)
})

test_that("reference: band is two-sided ([ref-crit, ref+crit])", {
  res <- make_compare_single(n = 10L, seed = 31L)
  ref_var <- as.character(res$Variety[5L])
  df      <- plot_compare(res, type = "dotplot", reference = ref_var,
                          return_data = TRUE)
  crit    <- unique(df$crit)
  expect_equal(unique(df$band_lo), unique(df$ref_pred) - crit)
  expect_equal(unique(df$band_hi), unique(df$ref_pred) + crit)
})

test_that("reference: sig has three levels (better/ns/worse) when spread is wide", {
  # Use a tiny criterion so varieties above AND below reference are significant
  res        <- make_compare_single(n = 20L, seed = 32L)
  res$HSD    <- rep(0.01, 20L)   # tiny criterion
  ref_var    <- as.character(res$Variety[10L])   # middle variety
  df         <- plot_compare(res, type = "dotplot", reference = ref_var,
                             return_data = TRUE)
  sig_vals   <- unique(df$sig)
  expect_true("better" %in% sig_vals)
  expect_true("worse"  %in% sig_vals)
})

test_that("reference: reference variety itself gets sig='ns'", {
  res     <- make_compare_single(n = 10L, seed = 33L)
  ref_var <- as.character(res$Variety[5L])
  df      <- plot_compare(res, type = "dotplot", reference = ref_var,
                          return_data = TRUE)
  expect_equal(df$sig[df$Variety == ref_var], "ns")
})

test_that("reference: unknown variety name warns and falls back to top-ranked", {
  res <- make_compare_single(n = 10L, seed = 34L)
  expect_warning(
    df <- plot_compare(res, type = "dotplot", reference = "BadVar",
                       return_data = TRUE),
    "not found in group"
  )
  # After fallback: ref_pred should equal max pred
  expect_equal(unique(df$ref_pred), max(df$pred))
})

test_that("reference: ignored with warning for type='letters'", {
  res <- make_compare_single(n = 8L, seed = 35L)
  expect_warning(
    p <- plot_compare(res, type = "letters", reference = "G01"),
    "only used for type"
  )
  expect_s3_class(p, "ggplot")
})

test_that("reference: ignored with warning for type='heatmap'", {
  res <- make_compare_single(n = 6L, seed = 36L)
  expect_warning(
    p <- plot_compare(res, type = "heatmap", reference = "G01"),
    "only used for type"
  )
  expect_s3_class(p, "ggplot")
})

test_that("reference: non-character stops", {
  res <- make_compare_single()
  expect_error(plot_compare(res, reference = 123L), "'reference' must be")
})

# ---------------------------------------------------------------------------
# 10. interactive argument
# ---------------------------------------------------------------------------

test_that("interactive=FALSE returns a ggplot (plotly not needed)", {
  res <- make_compare_single()
  p   <- plot_compare(res, interactive = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("interactive=TRUE returns a pc_interactive object", {
  res <- make_compare_single()
  p   <- plot_compare(res, interactive = TRUE)
  expect_s3_class(p, "pc_interactive")
})

test_that("pc_add() adds ggplot2 layers to a pc_interactive object", {
  res <- make_compare_single()
  p   <- pc_add(plot_compare(res, interactive = TRUE),
                ggplot2::ggtitle("Interactive title"))
  expect_s3_class(p, "pc_interactive")
  expect_equal(biomAid:::.pc_resolve_plot(p)$labels$title, "Interactive title")
})

test_that("interactive=TRUE errors at print if plotly not available", {
  # Skip if plotly is installed — we can only test the error path when absent
  skip_if(requireNamespace("plotly", quietly = TRUE),
          "plotly is installed; cannot test absence error")
  res <- make_compare_single()
  p   <- plot_compare(res, interactive = TRUE)
  expect_s3_class(p, "pc_interactive")        # construction succeeds without plotly
  expect_error(print(p), "plotly.*required")  # error fires at print time
})

test_that("interactive: non-logical stops", {
  res <- make_compare_single()
  expect_error(plot_compare(res, interactive = "yes"), "single logical")
})

# ---------------------------------------------------------------------------
# 11. LSD and Bonferroni criterion columns
# ---------------------------------------------------------------------------

test_that("plot_compare works with LSD criterion column", {
  res <- make_compare_single(type = "LSD", seed = 20L)
  expect_s3_class(plot_compare(res, type = "dotplot"), "ggplot")
})

test_that("plot_compare works with Bonferroni criterion column", {
  res <- make_compare_single(type = "Bonferroni", seed = 21L)
  expect_s3_class(plot_compare(res, type = "letters"), "ggplot")
})

# ---------------------------------------------------------------------------
# 10. ggplot extensibility
# ---------------------------------------------------------------------------

test_that("returned ggplot can be extended with +", {
  res <- make_compare_single()
  p   <- plot_compare(res) +
           ggplot2::ggtitle("Test title")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Test title")
})
