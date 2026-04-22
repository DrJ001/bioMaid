# tests/testthat/test-plot_randomRegress.R
# Tests for plot_randomRegress() — no ASReml needed.

library(ggplot2)

# ---- Shared mock result ------------------------------------------------
make_mock_res <- function(ns = 3L, nvar = 8L,
                          levs = c("T0", "T1", "T2")) {
  set.seed(42)
  sites <- paste0("E", seq_len(ns))

  blups <- do.call(rbind, lapply(sites, function(s) {
    u0 <- rnorm(nvar, 0, 200)
    u1 <- 0.9  * u0 + rnorm(nvar, 0, 80)
    u2 <- 1.2  * u0 + rnorm(nvar, 0, 120)
    data.frame(
      Site    = s,
      Variety = paste0("V", seq_len(nvar)),
      T0 = u0, T1 = u1, T2 = u2,
      resp.T1 = u1 - 0.9  * u0,
      resp.T2 = u2 - 1.2  * u0,
      HSD.T1  = NA_real_,
      HSD.T2  = NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  beta <- list(
    T1 = matrix(runif(ns, 0.7, 1.1), ns, 1L,
                dimnames = list(sites, "T0")),
    T2 = matrix(runif(ns, 1.0, 1.4), ns, 1L,
                dimnames = list(sites, "T0"))
  )

  nts  <- length(levs) * ns
  nms  <- paste0(rep(levs, ns), "-", rep(sites, each = length(levs)))
  G    <- crossprod(matrix(rnorm(nts * nts), nts, nts))
  dimnames(G) <- list(nms, nms)

  list(
    blups     = blups,
    beta      = beta,
    Gmat      = G,
    TGmat     = G,
    sigmat    = matrix(runif(ns * 2L, 100, 400), ns, 2L,
                       dimnames = list(sites, c("T1", "T2"))),
    tmat      = diag(nts),
    cond_list = list(T0 = NULL, T1 = "T0", T2 = "T0"),
    type      = "baseline"
  )
}

# ---- Input validation --------------------------------------------------

test_that("invalid type gives informative error", {
  expect_error(plot_randomRegress(make_mock_res(), type = "circles"),
               regexp = "arg")
})

test_that("missing res elements give informative error", {
  bad <- make_mock_res()
  bad$Gmat <- NULL
  expect_error(plot_randomRegress(bad), regexp = "Gmat")
})

test_that("non-theme theme argument gives informative error", {
  expect_error(
    plot_randomRegress(make_mock_res(), theme = "bw"),
    regexp = "theme"
  )
})

test_that("treatments with no match gives informative error", {
  expect_error(
    plot_randomRegress(make_mock_res(), type = "regress",
                       treatments = "N99"),
    regexp = "No conditioned treatments"
  )
})

# ---- Return types -------------------------------------------------------

test_that("all three types return a ggplot object", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant", "gmat"))
    expect_s3_class(plot_randomRegress(res, type = tp), "ggplot")
})

test_that("highlight = NULL suppresses highlighting without error", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant"))
    expect_s3_class(
      plot_randomRegress(res, type = tp, highlight = NULL), "ggplot"
    )
})

test_that("user-specified highlight returns ggplot", {
  res  <- make_mock_res()
  vars <- unique(res$blups$Variety)[1:3]
  expect_s3_class(
    plot_randomRegress(res, type = "quadrant", highlight = vars),
    "ggplot"
  )
})

test_that("non-existent highlight variety raises a warning", {
  res <- make_mock_res()
  expect_warning(
    plot_randomRegress(res, type = "quadrant",
                       highlight = c("V1", "DoesNotExist")),
    regexp = "not found"
  )
})

test_that("default highlights return data frame with Variety and group cols", {
  res   <- make_mock_res(ns = 4L, nvar = 12L)
  qdata <- plot_randomRegress(res, type = "quadrant", return_data = TRUE)
  hl    <- biomAid:::.rreg_default_highlights(qdata)
  expect_s3_class(hl, "data.frame")
  expect_true(all(c("Variety", "group") %in% names(hl)))
  expect_true(all(hl$group %in% c("tr", "bl")))
})

test_that("default highlights contain at most 6 varieties", {
  res   <- make_mock_res(ns = 4L, nvar = 12L)
  qdata <- plot_randomRegress(res, type = "quadrant", return_data = TRUE)
  hl    <- biomAid:::.rreg_default_highlights(qdata)
  expect_lte(nrow(hl), 6L)
  expect_gte(nrow(hl), 1L)
})

test_that("return_data = TRUE returns a data.frame for all types", {
  res <- make_mock_res()
  for (tp in c("regress", "quadrant", "gmat")) {
    df <- plot_randomRegress(res, type = tp, return_data = TRUE)
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0L)
  }
})

# ---- Data frame column checks ------------------------------------------

test_that("regress return_data has required columns", {
  df <- plot_randomRegress(make_mock_res(), type = "regress",
                           return_data = TRUE)
  expect_true(all(c("Site", "Variety", "x", "y",
                    "pair_label", "beta") %in% names(df)))
})

test_that("quadrant return_data has required columns", {
  df <- plot_randomRegress(make_mock_res(), type = "quadrant",
                           return_data = TRUE)
  expect_true(all(c("Site", "Variety", "x", "y",
                    "pair_label") %in% names(df)))
})

test_that("gmat return_data has row_var, col_var, corr", {
  df <- plot_randomRegress(make_mock_res(), type = "gmat",
                           return_data = TRUE)
  expect_true(all(c("row_var", "col_var", "corr") %in% names(df)))
})

# ---- Correctness checks ------------------------------------------------

test_that("gmat correlations are within [-1, 1]", {
  df <- plot_randomRegress(make_mock_res(), type = "gmat",
                           return_data = TRUE)
  expect_true(all(df$corr >= -1 - 1e-9 & df$corr <= 1 + 1e-9))
})

test_that("gmat diagonal correlations equal 1", {
  df <- plot_randomRegress(make_mock_res(), type = "gmat",
                           return_data = TRUE)
  diag_rows <- as.character(df$row_var) == as.character(df$col_var)
  expect_equal(df$corr[diag_rows], rep(1, sum(diag_rows)), tolerance = 1e-9)
})

test_that("regress data has one row per variety x site x conditioned treatment", {
  res <- make_mock_res(ns = 3L, nvar = 8L)
  df  <- plot_randomRegress(res, type = "regress", return_data = TRUE)
  # 2 conditioned treatments x 3 sites x 8 varieties = 48
  expect_equal(nrow(df), 2L * 3L * 8L)
})

test_that("treatments filter reduces rows in regress output", {
  res  <- make_mock_res(ns = 3L, nvar = 8L)
  df_all <- plot_randomRegress(res, type = "regress",
                               return_data = TRUE)
  df_sub <- plot_randomRegress(res, type = "regress",
                               treatments  = "T1",
                               return_data = TRUE)
  expect_lt(nrow(df_sub), nrow(df_all))
  expect_true(all(df_sub$pair_label == "T1 | T0"))
})

test_that("regress return_data contains beta column with finite values", {
  df <- plot_randomRegress(make_mock_res(), type = "regress",
                           return_data = TRUE)
  expect_true("beta" %in% names(df))
  expect_true(all(is.finite(df$beta)))
})

test_that("custom theme is applied without error", {
  res <- make_mock_res()
  expect_s3_class(
    plot_randomRegress(res, type = "quadrant",
                       theme = ggplot2::theme_classic()),
    "ggplot"
  )
})
