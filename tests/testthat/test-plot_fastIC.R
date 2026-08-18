# tests/testthat/test-plot_fastIC.R
# Tests for plot_fastIC() — the ggplot2 visualisation layer for fastIC().
#
# Structure
# =========
#  A. Input validation           – bad res, bad type, bad arguments.
#  B. Removed type rejection     – old names ("cve","trajectory","pairs",
#                                  "loads","specialist","winner") must error.
#  C. Precondition checks        – missing columns trigger informative errors.
#  D. Return value structure     – each type returns a ggplot / list correctly.
#  E. return_data = TRUE         – data component has expected columns/dims.
#  F. Highlight argument         – NULL, named vector, "default" all work.
#  G. Internal helpers           – .pfi_parse(), .pfi_geno_data(), etc.
#
# All tests use a synthetic fastIC() result built directly (no ASReml).

# ===========================================================================
# SHARED FIXTURE
# ===========================================================================
# Build a minimal but complete fastIC() result: 4 envs x 6 genotypes x 3 factors
# ic.num = 2 -> 2 iClasses ("pp" and "pn")

local({

  loads <- matrix(
    c( 2.0,  1.0,  0.4,
       1.5, -0.8,  0.6,
       1.8,  0.6, -0.5,
       1.2, -1.2,  0.3),
    nrow = 4L, ncol = 3L, byrow = TRUE,
    dimnames = list(c("E1","E2","E3","E4"), paste0("loads", 1:3))
  )

  scores <- matrix(
    c( 1.0,  0.5,  0.3,
       0.2, -0.3,  0.4,
      -0.5,  0.8, -0.2,
       0.7,  0.1,  0.6,
      -0.3, -0.6,  0.5,
       0.9, -0.2, -0.4),
    nrow = 6L, ncol = 3L, byrow = TRUE,
    dimnames = list(paste0("G", 1:6), paste0("score", 1:3))
  )

  spec_var <- c(E1 = 0.10, E2 = 0.20, E3 = 0.15, E4 = 0.25)
  envs     <- c("E1","E2","E3","E4")
  genos    <- paste0("G", 1:6)
  m        <- 6L; t <- 4L; k <- 3L

  CVE_mat    <- scores %*% t(loads)
  env_rep    <- rep(seq_len(t), each = m)
  geno_rep   <- rep(seq_len(m), times = t)

  df <- data.frame(
    Site     = factor(envs[env_rep],  levels = envs),
    Genotype = factor(genos[geno_rep], levels = genos),
    loads1   = loads[env_rep, 1L],
    loads2   = loads[env_rep, 2L],
    loads3   = loads[env_rep, 3L],
    spec.var = spec_var[env_rep],
    score1   = scores[geno_rep, 1L],
    score2   = scores[geno_rep, 2L],
    score3   = scores[geno_rep, 3L],
    CVE      = CVE_mat[cbind(geno_rep, env_rep)],
    stringsAsFactors = FALSE
  )

  for (r in 1:k)
    df[[paste0("fitted", r)]] <- df[[paste0("score", r)]] *
                                  df[[paste0("loads", r)]]

  # global_op / global_dev / global_stab
  fitted1_mat  <- outer(scores[, 1L], loads[, 1L])
  dev_mat      <- CVE_mat - fitted1_mat
  df$global_op   <- mean(loads[, 1L]) * df$score1
  df$global_dev  <- dev_mat[cbind(geno_rep, env_rep)]
  df$global_stab <- sqrt(rowMeans(dev_mat^2))[geno_rep]

  # iClass (ic.num = 2)
  sign_str <- apply(loads[, 1:2], 1L, function(x)
    paste(ifelse(x >= 0, "p", "n"), collapse = ""))
  ic_levs  <- unique(sign_str)
  df$iclass <- factor(sign_str[env_rep], levels = ic_levs)

  # iClassOP and iClassRMSD
  fitted_ic <- outer(scores[, 1L], loads[, 1L]) +
               outer(scores[, 2L], loads[, 2L])
  iClassOP_mat   <- matrix(NA_real_, m, length(ic_levs),
                            dimnames = list(genos, ic_levs))
  iClassRMSD_mat <- matrix(NA_real_, m, length(ic_levs),
                            dimnames = list(genos, ic_levs))
  for (w in ic_levs) {
    env_w  <- which(sign_str == w)
    mld_w  <- colMeans(loads[env_w, 1:2, drop = FALSE])
    iClassOP_mat[, w]   <- scores[, 1:2, drop = FALSE] %*% mld_w
    dev_w               <- (CVE_mat - fitted_ic)[, env_w, drop = FALSE]
    iClassRMSD_mat[, w] <- sqrt(rowMeans(dev_w^2))
  }
  gi <- cbind(as.character(df$Genotype), as.character(df$iclass))
  df$iClassOP   <- iClassOP_mat[gi]
  df$iClassRMSD <- iClassRMSD_mat[gi]

  df <- df[order(df$iclass, df$Site, df$Genotype), ]
  rownames(df) <- NULL

  # ---- VAF attributes (mirrors fastIC() computation) --------------------
  loads_sq  <- loads^2                              # t x k
  total_var <- rowSums(loads_sq) + spec_var         # length t
  vaf_mat   <- loads_sq / total_var                 # t x k proportions

  vaf_env_df        <- as.data.frame(vaf_mat)
  names(vaf_env_df) <- paste0("Factor", 1:k)
  vaf_env_df$Site      <- envs
  vaf_env_df$Specific  <- spec_var / total_var
  vaf_env_df$total_var <- total_var
  vaf_env_df <- vaf_env_df[, c("Site", paste0("Factor", 1:k), "Specific", "total_var")]
  rownames(vaf_env_df) <- NULL

  total_all  <- sum(total_var)
  pct_factor <- colSums(loads_sq) / total_all
  pct_spec   <- sum(spec_var) / total_all
  vaf_summ   <- data.frame(
    factor  = c(paste0("Factor ", 1:k), "Specific"),
    pct_var = c(pct_factor, pct_spec),
    stringsAsFactors = FALSE
  )
  vaf_summ$cum_pct <- cumsum(vaf_summ$pct_var)

  attr(df, "vaf_env")     <- vaf_env_df
  attr(df, "vaf_summary") <- vaf_summ

  # Assign to parent env so tests can use .pfi_res
  .pfi_res <<- df
})

# ===========================================================================
# SECTION A – Input validation
# ===========================================================================

test_that("A: non-data-frame res errors", {
  expect_error(plot_fastIC("not a df"), "data.frame")
})

test_that("A: res with < 4 columns errors", {
  expect_error(plot_fastIC(data.frame(a = 1, b = 2, c = 3)), "data.frame")
})

test_that("A: invalid type errors with informative message", {
  expect_error(plot_fastIC(.pfi_res, type = "banana"), "should be one of")
})

test_that("A: non-theme theme argument errors", {
  expect_error(
    plot_fastIC(.pfi_res, type = "fast", theme = "not_a_theme"),
    "theme"
  )
})

test_that("A: non-logical return_data errors", {
  expect_error(
    plot_fastIC(.pfi_res, type = "fast", return_data = "yes"),
    "return_data"
  )
})

# ===========================================================================
# SECTION B – Removed type names are rejected
# ===========================================================================

test_that("B: old type 'cve' is rejected (renamed to 'CVE')", {
  expect_error(plot_fastIC(.pfi_res, type = "cve"), "should be one of")
})

test_that("B: old type 'trajectory' is rejected (renamed to 'OP.variety')", {
  expect_error(plot_fastIC(.pfi_res, type = "trajectory"), "should be one of")
})

test_that("B: old type 'pairs' is rejected (renamed to 'OP.pairs')", {
  expect_error(plot_fastIC(.pfi_res, type = "pairs"), "should be one of")
})

test_that("B: removed type 'loads' is rejected", {
  expect_error(plot_fastIC(.pfi_res, type = "loads"), "should be one of")
})

test_that("B: removed type 'specialist' is rejected", {
  expect_error(plot_fastIC(.pfi_res, type = "specialist"), "should be one of")
})

test_that("B: removed type 'winner' is rejected", {
  expect_error(plot_fastIC(.pfi_res, type = "winner"), "should be one of")
})

test_that("B: old type 'dev' is rejected (removed)", {
  expect_error(plot_fastIC(.pfi_res, type = "dev"), "should be one of")
})

test_that("B: old type 'VE' is rejected (replaced by 'VAF')", {
  expect_error(plot_fastIC(.pfi_res, type = "VE"), "should be one of")
})

test_that("B: 'VAF' is rejected (moved to plot_faSummary())", {
  expect_error(plot_fastIC(.pfi_res, type = "VAF"), "should be one of")
})

# ===========================================================================
# SECTION C – Precondition checks
# ===========================================================================

# Strip a column from a copy of the fixture to trigger precondition errors

test_that("C: 'fast' errors when global_op is absent", {
  res_no_op <- .pfi_res[, setdiff(names(.pfi_res), "global_op")]
  expect_error(plot_fastIC(res_no_op, type = "fast"), "global_op")
})

test_that("C: 'fast' errors when global_stab is absent", {
  res_no_stab <- .pfi_res[, setdiff(names(.pfi_res), "global_stab")]
  expect_error(plot_fastIC(res_no_stab, type = "fast"), "global_stab")
})

test_that("C: 'iclass' errors when iclass column is absent", {
  res_no_ic <- .pfi_res[, setdiff(names(.pfi_res), c("iclass","iClassOP","iClassRMSD"))]
  expect_error(plot_fastIC(res_no_ic, type = "iclass"), "iclass")
})

test_that("C: 'OP.variety' errors when iclass column is absent", {
  res_no_ic <- .pfi_res[, setdiff(names(.pfi_res), c("iclass","iClassOP","iClassRMSD"))]
  expect_error(plot_fastIC(res_no_ic, type = "OP.variety"), "iclass")
})

test_that("C: 'OP.pairs' errors when iclass column is absent", {
  res_no_ic <- .pfi_res[, setdiff(names(.pfi_res), c("iclass","iClassOP","iClassRMSD"))]
  expect_error(plot_fastIC(res_no_ic, type = "OP.pairs"), "iclass")
})

test_that("C: 'biplot' errors when k < 2", {
  # Build a k=1 fixture
  res_1f <- .pfi_res[, setdiff(names(.pfi_res), c("loads2","loads3","score2","score3","fitted2","fitted3"))]
  expect_error(plot_fastIC(res_1f, type = "biplot"), "k >= 2")
})

test_that("C: biplot_factors wrong length errors", {
  expect_error(
    plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(1L, 2L, 3L)),
    "length-2"
  )
})

test_that("C: biplot_factors duplicate values errors", {
  expect_error(
    plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(2L, 2L)),
    "distinct"
  )
})

test_that("C: biplot_factors out-of-range errors", {
  expect_error(
    plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(1L, 5L)),
    "k = 3"
  )
})

test_that("C: biplot_factors = c(1,3) produces a ggplot (k=3 fixture)", {
  p <- plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(1L, 3L))
  expect_s3_class(p, "ggplot")
})

test_that("C: biplot_factors axis labels reflect chosen factors", {
  p <- plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(1L, 3L))
  expect_match(p$labels$x, "Factor 1")
  expect_match(p$labels$y, "Factor 3")
})

test_that("C: biplot_factors title notes non-default axes", {
  p <- plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(2L, 3L))
  expect_match(p$labels$title, "Factors 2 & 3")
})

test_that("C: biplot_factors default c(1,2) title has no factor annotation", {
  p <- plot_fastIC(.pfi_res, type = "biplot", biplot_factors = c(1L, 2L))
  expect_false(grepl("Factors [0-9]", p$labels$title))
})

test_that("C: biplot_factors ignored for non-biplot types", {
  # Should not error even with an out-of-range value when type != biplot
  expect_s3_class(
    plot_fastIC(.pfi_res, type = "fast", biplot_factors = c(9L, 99L)),
    "ggplot"
  )
})

# ===========================================================================
# SECTION D – Return value: ggplot object or list
# ===========================================================================

test_that("D: 'fast' returns a ggplot when return_data = FALSE", {
  p <- plot_fastIC(.pfi_res, type = "fast")
  expect_s3_class(p, "ggplot")
})

test_that("D: 'biplot' returns a ggplot", {
  p <- plot_fastIC(.pfi_res, type = "biplot")
  expect_s3_class(p, "ggplot")
})

test_that("D: 'CVE' returns a ggplot", {
  p <- plot_fastIC(.pfi_res, type = "CVE")
  expect_s3_class(p, "ggplot")
})

test_that("D: 'iclass' returns a ggplot", {
  p <- plot_fastIC(.pfi_res, type = "iclass")
  expect_s3_class(p, "ggplot")
})

test_that("D: 'OP.variety' returns a ggplot", {
  p <- plot_fastIC(.pfi_res, type = "OP.variety")
  expect_s3_class(p, "ggplot")
})

test_that("D: 'OP.pairs' returns a ggplot", {
  p <- plot_fastIC(.pfi_res, type = "OP.pairs")
  expect_s3_class(p, "ggplot")
})

test_that("D: return_data = TRUE returns a named list with $plot and $data", {
  for (typ in c("fast", "biplot", "CVE", "iclass", "OP.variety", "OP.pairs")) {
    out <- plot_fastIC(.pfi_res, type = typ, return_data = TRUE)
    expect_type(out, "list")
    expect_true("plot" %in% names(out), label = paste(typ, "has $plot"))
    expect_true("data" %in% names(out), label = paste(typ, "has $data"))
    expect_true(inherits(out$plot, "ggplot"), label = paste(typ, "$plot is ggplot"))
  }
})

test_that("D: default type is 'fast'", {
  p <- plot_fastIC(.pfi_res)
  expect_s3_class(p, "ggplot")
  # Verify it is the FAST plot by checking axis labels
  built <- ggplot2::ggplot_build(p)
  expect_true(
    grepl("global_op", built$plot$labels$x, ignore.case = TRUE) ||
    grepl("Overall", built$plot$labels$x, ignore.case = TRUE)
  )
})

# ===========================================================================
# SECTION E – return_data = TRUE: data component structure
# ===========================================================================

test_that("E: 'fast' return_data has global_op and global_stab columns", {
  out <- plot_fastIC(.pfi_res, type = "fast", return_data = TRUE)
  expect_true("global_op"   %in% names(out$data))
  expect_true("global_stab" %in% names(out$data))
  expect_true("highlighted"  %in% names(out$data))
})

test_that("E: 'fast' return_data has one row per genotype", {
  out <- plot_fastIC(.pfi_res, type = "fast", return_data = TRUE)
  n_geno <- length(unique(.pfi_res$Genotype))
  expect_equal(nrow(out$data), n_geno)
})

test_that("E: 'CVE' return_data is a wide matrix (genotype rows, env columns)", {
  out <- plot_fastIC(.pfi_res, type = "CVE", return_data = TRUE)
  expect_true(is.data.frame(out$data))
  # Should have genotype column + one column per environment
  n_envs <- length(unique(.pfi_res$Site))
  expect_equal(ncol(out$data), n_envs + 1L)
})

test_that("E: 'iclass' return_data has iclass, iClassOP, iClassRMSD columns", {
  out <- plot_fastIC(.pfi_res, type = "iclass", return_data = TRUE)
  expect_true("iclass"     %in% names(out$data))
  expect_true("iClassOP"   %in% names(out$data))
  expect_true("iClassRMSD" %in% names(out$data))
})

test_that("E: 'OP.variety' return_data has iclass and iClassOP columns", {
  out <- plot_fastIC(.pfi_res, type = "OP.variety", return_data = TRUE)
  expect_true("iclass"   %in% names(out$data))
  expect_true("iClassOP" %in% names(out$data))
})

test_that("E: 'OP.pairs' return_data has x_val, y_val, x_class, y_class", {
  out <- plot_fastIC(.pfi_res, type = "OP.pairs", return_data = TRUE)
  for (col in c("x_val", "y_val", "x_class", "y_class"))
    expect_true(col %in% names(out$data),
                label = paste("OP.pairs data has column", col))
})

test_that("E: 'OP.pairs' return_data is lower-triangular (i > j pairs only)", {
  out <- plot_fastIC(.pfi_res, type = "OP.pairs", return_data = TRUE)
  # With 2 iClasses there is exactly 1 pair
  n_ic    <- length(levels(.pfi_res$iclass))
  n_pairs <- n_ic * (n_ic - 1L) / 2L
  unique_pairs <- unique(out$data[, c("x_class","y_class")])
  expect_equal(nrow(unique_pairs), n_pairs)
})

test_that("E: 'biplot' return_data is a list with $geno and $envs", {
  out <- plot_fastIC(.pfi_res, type = "biplot", return_data = TRUE)
  expect_type(out$data, "list")
  expect_true("geno" %in% names(out$data))
  expect_true("envs" %in% names(out$data))
})

# ===========================================================================
# SECTION F – Highlight argument
# ===========================================================================

test_that("F: highlight = NULL produces a ggplot without error", {
  p <- plot_fastIC(.pfi_res, type = "fast", highlight = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("F: highlight = named character vector is accepted", {
  p <- plot_fastIC(.pfi_res, type = "fast", highlight = c("G1", "G3"))
  expect_s3_class(p, "ggplot")
})

test_that("F: highlight = 'default' produces a ggplot without error", {
  for (typ in c("fast", "biplot", "CVE", "iclass", "OP.variety", "OP.pairs")) {
    p <- plot_fastIC(.pfi_res, type = typ, highlight = "default")
    expect_true(inherits(p, "ggplot"), label = paste(typ, "default highlight is ggplot"))
  }
})

test_that("F: n_highlight = 1 is accepted", {
  p <- plot_fastIC(.pfi_res, type = "fast", n_highlight = 1L)
  expect_s3_class(p, "ggplot")
})

test_that("F: highlight names not in data are silently ignored", {
  # Only "G1" is valid; "ZZZZ" is not — should not error
  p <- plot_fastIC(.pfi_res, type = "fast", highlight = c("G1", "ZZZZ"))
  expect_s3_class(p, "ggplot")
})

# ===========================================================================
# SECTION G – Internal helper unit tests
# ===========================================================================

test_that("G: .pfi_parse() identifies sterm and gterm correctly", {
  p <- biomAid:::.pfi_parse(.pfi_res)
  expect_equal(p$sterm, "Site")
  expect_equal(p$gterm, "Genotype")
})

test_that("G: .pfi_parse() counts k = 3 factors", {
  p <- biomAid:::.pfi_parse(.pfi_res)
  expect_equal(p$k, 3L)
})

test_that("G: .pfi_parse() sets has_fast, has_stab, has_iclass TRUE", {
  p <- biomAid:::.pfi_parse(.pfi_res)
  expect_true(p$has_fast)
  expect_true(p$has_stab)
  expect_true(p$has_iclass)
})

test_that("G: .pfi_parse() sets has_fast FALSE when global_op absent", {
  res2 <- .pfi_res[, setdiff(names(.pfi_res), "global_op")]
  p    <- biomAid:::.pfi_parse(res2)
  expect_false(p$has_fast)
})

test_that("G: 'vaf_env' attribute is still carried on fastIC() output", {
  # Ownership of the VAF calculation moved to faSummary(); fastIC() passes it
  # through so plot_faSummary(type = "VAF") can consume it.
  expect_false(is.null(attr(.pfi_res, "vaf_env")))
})

test_that("G: 'vaf_summary' has 'factor', 'pct_var', 'cum_pct' columns", {
  summ <- attr(.pfi_res, "vaf_summary")
  expect_true("factor"  %in% names(summ))
  expect_true("pct_var" %in% names(summ))
  expect_true("cum_pct" %in% names(summ))
})

test_that("G: .pfi_geno_data() returns one row per genotype", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  expect_equal(nrow(gd), length(unique(.pfi_res$Genotype)))
})

test_that("G: .pfi_geno_data() contains global_op and global_stab", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  expect_true("global_op"   %in% names(gd))
  expect_true("global_stab" %in% names(gd))
})

test_that("G: .pfi_geno_data() pivots iClassOP wide (one col per iClass)", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  ic_levels <- levels(.pfi_res$iclass)
  for (ic in ic_levels)
    expect_true(paste0("iClassOP_", ic) %in% names(gd),
                label = paste("iClassOP_", ic, "present"))
})

test_that("G: .pfi_env_data() returns one row per environment", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  ed <- biomAid:::.pfi_env_data(.pfi_res, p)
  expect_equal(nrow(ed), length(unique(.pfi_res$Site)))
})

test_that("G: .pfi_highlights() returns character(0) for highlight = NULL", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  hl <- biomAid:::.pfi_highlights(gd, p$gterm, NULL, 3L, "fast")
  expect_equal(hl, character(0L))
})

test_that("G: .pfi_highlights() returns named varieties for explicit vector", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  hl <- biomAid:::.pfi_highlights(gd, p$gterm, c("G1","G2"), 3L, "fast")
  expect_equal(sort(hl), c("G1","G2"))
})

test_that("G: .pfi_highlights() default for 'fast' returns <= 2*n_highlight varieties", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  hl <- biomAid:::.pfi_highlights(gd, p$gterm, "default", 2L, "fast")
  expect_lte(length(hl), 4L)   # up to 2 * n_highlight
  expect_gte(length(hl), 1L)
})

test_that("G: .pfi_highlights() default for 'OP.variety' returns <= n_highlight varieties", {
  p  <- biomAid:::.pfi_parse(.pfi_res)
  gd <- biomAid:::.pfi_geno_data(.pfi_res, p)
  hl <- biomAid:::.pfi_highlights(gd, p$gterm, "default", 3L, "OP.variety")
  expect_lte(length(hl), 3L)
  expect_gte(length(hl), 1L)
})
