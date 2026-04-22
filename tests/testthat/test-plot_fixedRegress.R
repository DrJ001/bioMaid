# tests/testthat/test-plot_fixedRegress.R
# Tests for plot_fixedRegress() -- no ASReml needed.

library(ggplot2)

# ---- Shared mock result ------------------------------------------------
make_mock_freg <- function(ns = 4L, nvar = 12L,
                            levs = c("T0", "T1", "T2"),
                            by_col   = "Site",
                            geno_col = "Genotype") {
  set.seed(7)
  groups <- paste0("G", seq_len(ns))

  blues_list <- lapply(groups, function(g) {
    t0 <- rnorm(nvar, 0, 200)
    t1 <- 0.9  * t0 + 50 + rnorm(nvar, 0, 80)
    t2 <- 1.2  * t0 + 100 + rnorm(nvar, 0, 120)
    r1 <- residuals(lm(t1 ~ t0))
    r2 <- residuals(lm(t2 ~ t0))
    df <- data.frame(t0, t1, t2, r1, r2,
                     se1 = abs(rnorm(nvar, 10, 2)),
                     se2 = abs(rnorm(nvar, 12, 2)),
                     HSD1 = abs(rnorm(1, 50, 5)),
                     HSD2 = abs(rnorm(1, 55, 5)),
                     stringsAsFactors = FALSE)
    names(df) <- c(levs, paste0("resp.", levs[-1]),
                   paste0("se.", levs[-1]), paste0("HSD.", levs[-1]))
    cbind(
      setNames(data.frame(g, paste0("V", sprintf("%02d", seq_len(nvar))),
                          stringsAsFactors = FALSE),
               c(by_col, geno_col)),
      df
    )
  })
  blues <- do.call(rbind, blues_list)

  # beta: named list, one data frame per conditioned treatment
  beta <- list(
    T1 = data.frame(
      G = groups,
      T0 = runif(ns, 0.7, 1.1),
      stringsAsFactors = FALSE
    ),
    T2 = data.frame(
      G = groups,
      T0 = runif(ns, 1.0, 1.4),
      stringsAsFactors = FALSE
    )
  )
  names(beta$T1)[1] <- by_col
  names(beta$T2)[1] <- by_col

  sigmat <- matrix(runif(ns * 2L, 50, 150), ns, 2L,
                   dimnames = list(groups, c("T1", "T2")))

  list(
    blues     = blues,
    beta      = beta,
    sigmat    = sigmat,
    cond_list = list(T0 = NULL, T1 = "T0", T2 = "T0"),
    type      = "baseline"
  )
}

# ---- Input validation --------------------------------------------------

test_that("invalid type gives informative error", {
  expect_error(plot_fixedRegress(make_mock_freg(), type = "gmat"),
               regexp = "arg")
})

test_that("missing res elements give informative error", {
  bad <- make_mock_freg()
  bad$blues <- NULL
  expect_error(plot_fixedRegress(bad), regexp = "blues")
})

test_that("non-theme theme argument gives informative error", {
  expect_error(
    plot_fixedRegress(make_mock_freg(), theme = "bw"),
    regexp = "theme"
  )
})

test_that("treatments with no match gives informative error", {
  expect_error(
    plot_fixedRegress(make_mock_freg(), type = "regress",
                      treatments = "N99"),
    regexp = "No conditioned treatments"
  )
})

# ---- Return types -------------------------------------------------------

test_that("both types return a ggplot object", {
  res <- make_mock_freg()
  for (tp in c("regress", "quadrant"))
    expect_s3_class(plot_fixedRegress(res, type = tp), "ggplot")
})

test_that("return_data = TRUE returns a data.frame for both types", {
  res <- make_mock_freg()
  for (tp in c("regress", "quadrant")) {
    df <- plot_fixedRegress(res, type = tp, return_data = TRUE)
    expect_s3_class(df, "data.frame")
    expect_gt(nrow(df), 0L)
  }
})

# ---- Data frame column checks ------------------------------------------

test_that("regress return_data has required columns including intercept", {
  df <- plot_fixedRegress(make_mock_freg(), type = "regress",
                          return_data = TRUE)
  expect_true(all(c("Group", "Genotype", "x", "y",
                    "pair_label", "beta", "intercept") %in% names(df)))
})

test_that("quadrant return_data has required columns", {
  df <- plot_fixedRegress(make_mock_freg(), type = "quadrant",
                          return_data = TRUE)
  expect_true(all(c("Group", "Genotype", "x", "y",
                    "pair_label") %in% names(df)))
})

test_that("Group column contains the original group values", {
  res <- make_mock_freg(ns = 3L)
  df  <- plot_fixedRegress(res, type = "quadrant", return_data = TRUE)
  expect_true(all(c("G1", "G2", "G3") %in% unique(df$Group)))
})

test_that("intercept values are finite", {
  df <- plot_fixedRegress(make_mock_freg(), type = "regress",
                          return_data = TRUE)
  expect_true(all(is.finite(df$intercept)))
})

# ---- Highlight tests ---------------------------------------------------

test_that("highlight = NULL suppresses highlighting without error", {
  res <- make_mock_freg()
  for (tp in c("regress", "quadrant"))
    expect_s3_class(
      plot_fixedRegress(res, type = tp, highlight = NULL), "ggplot"
    )
})

test_that("user-specified highlight returns ggplot", {
  res   <- make_mock_freg()
  genos <- unique(res$blues$Genotype)[1:3]
  expect_s3_class(
    plot_fixedRegress(res, type = "quadrant", highlight = genos),
    "ggplot"
  )
})

test_that("non-existent highlight genotype raises a warning", {
  res <- make_mock_freg()
  expect_warning(
    plot_fixedRegress(res, type = "quadrant",
                      highlight = c("V01", "DoesNotExist")),
    regexp = "not found"
  )
})

test_that("default highlights return df with Genotype and group cols", {
  res   <- make_mock_freg(ns = 4L, nvar = 15L)
  qdata <- plot_fixedRegress(res, type = "quadrant", return_data = TRUE)
  hl    <- biomAid:::.freg_default_highlights(qdata)
  expect_s3_class(hl, "data.frame")
  expect_true(all(c("Genotype", "group") %in% names(hl)))
  expect_true(all(hl$group %in% c("tr", "bl")))
})

test_that("default highlights contain at most 6 genotypes", {
  res   <- make_mock_freg(ns = 4L, nvar = 15L)
  qdata <- plot_fixedRegress(res, type = "quadrant", return_data = TRUE)
  hl    <- biomAid:::.freg_default_highlights(qdata)
  expect_lte(nrow(hl), 6L)
  expect_gte(nrow(hl), 1L)
})

test_that("treatments filter reduces panels in regress output", {
  res    <- make_mock_freg()
  df_all <- plot_fixedRegress(res, type = "regress", return_data = TRUE)
  df_sub <- plot_fixedRegress(res, type = "regress",
                               treatments = "T1", return_data = TRUE)
  expect_lt(nrow(df_sub), nrow(df_all))
  expect_true(all(df_sub$pair_label == "T1 | T0"))
})

test_that("custom theme is applied without error", {
  res <- make_mock_freg()
  expect_s3_class(
    plot_fixedRegress(res, type = "quadrant",
                      theme = ggplot2::theme_classic()),
    "ggplot"
  )
})

test_that("plot works with by = NULL (single Group = 'All')", {
  res        <- make_mock_freg(ns = 1L, by_col = "Group")
  res$blues[["Group"]] <- "All"
  res$beta$T1[["Group"]] <- "All"
  res$beta$T2[["Group"]] <- "All"
  expect_s3_class(plot_fixedRegress(res, type = "regress"), "ggplot")
})
