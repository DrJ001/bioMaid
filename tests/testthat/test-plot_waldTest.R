# tests/testthat/test-plot_waldTest.R
# Tests for plot_waldTest() — pure R, no ASReml needed.

# ---------------------------------------------------------------------------
# Helper: build a waldTest result from a synthetic pred object.
# ---------------------------------------------------------------------------
make_wt_result <- function(treatments = c("N0", "N1", "N2"),
                            by_var     = NULL,
                            n_groups   = 2L,
                            include_zero = FALSE,
                            seed       = 42L) {
  set.seed(seed)
  n <- length(treatments)

  make_one_pred <- function(grp = NULL, s = seed) {
    set.seed(s)
    pv <- data.frame(
      Treatment       = factor(treatments, levels = treatments),
      predicted.value = stats::rnorm(n, 50, 10),
      std.error       = stats::runif(n, 1, 3),
      status          = factor(rep("Estimable", n)),
      stringsAsFactors = FALSE
    )
    if (!is.null(grp) && !is.null(by_var)) pv[[by_var]] <- factor(grp)
    A    <- matrix(stats::rnorm(n * n) * 0.3, n, n)
    vcov <- crossprod(A) + diag(stats::runif(n, 0.5, 1.5))
    list(pvals = pv, vcov = vcov)
  }

  if (!is.null(by_var)) {
    grp_names <- paste0("G", seq_len(n_groups))
    preds     <- mapply(make_one_pred, grp_names,
                        seq_len(n_groups), SIMPLIFY = FALSE)
    all_pv    <- do.call(rbind, lapply(preds, `[[`, "pvals"))
    n_tot     <- nrow(all_pv)
    vcov      <- matrix(0, n_tot, n_tot)
    for (i in seq_len(n_groups)) {
      idx <- ((i - 1L) * n + 1L):(i * n)
      vcov[idx, idx] <- preds[[i]]$vcov
    }
    pred <- list(pvals = all_pv, vcov = vcov)
  } else {
    pred <- make_one_pred()
  }

  cc <- list(list(coef = treatments, type = "con", comp = "pairwise"))
  if (include_zero)
    cc <- c(cc, list(list(coef  = treatments[2:3],
                          type  = "zero",
                          group = "joint_test")))

  waldTest(pred, cc = cc, by = by_var)
}

# Helper: Zero-only result (no contrasts)
make_zero_only <- function(seed = 42L) {
  set.seed(seed)
  treatments <- c("N0", "N1", "N2")
  n  <- 3L
  pv <- data.frame(
    Treatment       = factor(treatments),
    predicted.value = stats::rnorm(n, 50, 10),
    std.error       = stats::runif(n, 1, 3),
    status          = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(stats::rnorm(n * n) * 0.3, n, n)
  vcov <- crossprod(A) + diag(stats::runif(n, 0.5, 1.5))
  waldTest(list(pvals = pv, vcov = vcov),
           cc = list(list(coef  = c("N1", "N2"),
                          type  = "zero",
                          group = "N1:N2 joint")))
}

# ---------------------------------------------------------------------------
# 1. Basic return types
# ---------------------------------------------------------------------------
test_that("plot_waldTest returns a ggplot object", {
  res <- make_wt_result()
  p   <- plot_waldTest(res)
  expect_s3_class(p, "ggplot")
})

test_that("return_data = TRUE returns a data frame", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_s3_class(df, "data.frame")
})

test_that("plot_waldTest is returned invisibly", {
  res <- make_wt_result()
  expect_invisible(plot_waldTest(res))
})

# ---------------------------------------------------------------------------
# 2. Data frame structure (return_data = TRUE)
# ---------------------------------------------------------------------------
test_that("return_data has all required columns", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expected <- c("label", "group", "Estimate", "Std.Error",
                "CI_lower", "CI_upper", "x_label",
                "P.Value", "neg_log10_p", "significant",
                "p_label", "row_id", "y_factor")
  expect_true(all(expected %in% names(df)))
})

test_that("CI_lower < Estimate < CI_upper for all rows", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_true(all(df$CI_lower < df$Estimate))
  expect_true(all(df$Estimate < df$CI_upper))
})

test_that("neg_log10_p is non-negative for all rows", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_true(all(df$neg_log10_p >= 0))
})

test_that("significant column is logical", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_type(df$significant, "logical")
})

test_that("p_label starts with 'p<' or 'p='", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_true(all(grepl("^p[<=]", df$p_label)))
})

test_that("y_factor is a factor", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_s3_class(df$y_factor, "factor")
})

test_that("row_id values are unique", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_equal(length(unique(df$row_id)), nrow(df))
})

test_that("x_label equals CI_upper for all rows", {
  res <- make_wt_result()
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_equal(df$x_label, df$CI_upper)
})

test_that("nrow equals number of pairwise contrasts", {
  res <- make_wt_result(treatments = c("N0","N1","N2"))
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_equal(nrow(df), choose(3L, 2L))
})

# ---------------------------------------------------------------------------
# 3. ci_level
# ---------------------------------------------------------------------------
test_that("wider ci_level produces wider CI arms", {
  res  <- make_wt_result()
  df95 <- plot_waldTest(res, ci_level = 0.95, return_data = TRUE)
  df99 <- plot_waldTest(res, ci_level = 0.99, return_data = TRUE)
  w95  <- mean(df95$CI_upper - df95$CI_lower, na.rm = TRUE)
  w99  <- mean(df99$CI_upper - df99$CI_lower, na.rm = TRUE)
  expect_gt(w99, w95)
})

test_that("ci_level = 0.90 is accepted", {
  res <- make_wt_result()
  expect_no_error(plot_waldTest(res, ci_level = 0.90, return_data = TRUE))
})

test_that("ci_level = 0.99 is accepted", {
  res <- make_wt_result()
  expect_no_error(plot_waldTest(res, ci_level = 0.99, return_data = TRUE))
})

# ---------------------------------------------------------------------------
# 4. alpha
# ---------------------------------------------------------------------------
test_that("laxer alpha flags more rows as significant", {
  res  <- make_wt_result(seed = 99L)
  df1  <- plot_waldTest(res, alpha = 0.05, return_data = TRUE)
  df2  <- plot_waldTest(res, alpha = 0.50, return_data = TRUE)
  expect_gte(sum(df2$significant), sum(df1$significant))
})

test_that("alpha = 0.01 is accepted", {
  res <- make_wt_result()
  expect_no_error(plot_waldTest(res, alpha = 0.01, return_data = TRUE))
})

# ---------------------------------------------------------------------------
# 5. Zero-only result stops with informative error
# ---------------------------------------------------------------------------
test_that("Zero-only result raises a stop() error", {
  res <- make_zero_only()
  expect_error(plot_waldTest(res),
               "contains only joint zero test results")
})

test_that("Zero-only stop() message mentions 'con' contrasts", {
  res <- make_zero_only()
  expect_error(plot_waldTest(res), "'con' type contrasts")
})

# ---------------------------------------------------------------------------
# 6. Mixed Contrast + Zero result warns and plots contrasts only
# ---------------------------------------------------------------------------
test_that("mixed result raises a warning about Zero results", {
  res <- make_wt_result(include_zero = TRUE)
  expect_warning(plot_waldTest(res),
                 "zero test results.*cannot be displayed")
})

test_that("mixed result still returns a ggplot after warning", {
  res <- make_wt_result(include_zero = TRUE)
  p   <- suppressWarnings(plot_waldTest(res))
  expect_s3_class(p, "ggplot")
})

test_that("mixed result data contains only Contrast rows", {
  res <- make_wt_result(include_zero = TRUE)
  df  <- suppressWarnings(plot_waldTest(res, return_data = TRUE))
  # No Zero rows — only pairwise contrast rows
  expect_equal(nrow(df), nrow(res$Contrasts))
})

# ---------------------------------------------------------------------------
# 7. by-group / faceting
# ---------------------------------------------------------------------------
test_that("by-group with facet = TRUE returns ggplot", {
  res <- make_wt_result(by_var = "Site", n_groups = 3L)
  expect_s3_class(plot_waldTest(res, facet = TRUE), "ggplot")
})

test_that("by-group with facet = FALSE returns ggplot", {
  res <- make_wt_result(by_var = "Site", n_groups = 3L)
  expect_s3_class(plot_waldTest(res, facet = FALSE), "ggplot")
})

test_that("by-group data contains expected group labels", {
  res <- make_wt_result(by_var = "Site", n_groups = 2L)
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_true(all(c("G1", "G2") %in% df$group))
})

test_that("no-group result has group = 'All'", {
  res <- make_wt_result(by_var = NULL)
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_true(all(df$group == "All"))
})

test_that("row_ids are unique across groups", {
  res <- make_wt_result(by_var = "Site", n_groups = 2L)
  df  <- plot_waldTest(res, return_data = TRUE)
  expect_equal(length(unique(df$row_id)), nrow(df))
})

test_that("by-group nrow = n_groups * contrasts_per_group", {
  n_trt <- 3L
  n_grp <- 2L
  res   <- make_wt_result(by_var = "Site", n_groups = n_grp)
  df    <- plot_waldTest(res, return_data = TRUE)
  expect_equal(nrow(df), choose(n_trt, 2L) * n_grp)
})

# ---------------------------------------------------------------------------
# 8. F-test result
# ---------------------------------------------------------------------------
test_that("F-test waldTest result plots correctly", {
  set.seed(7L)
  treatments <- c("N0", "N1", "N2")
  n  <- 3L
  pv <- data.frame(
    Treatment       = factor(treatments),
    predicted.value = stats::rnorm(n, 50, 10),
    std.error       = stats::runif(n, 1, 3),
    status          = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(stats::rnorm(n * n) * 0.3, n, n)
  vcov <- crossprod(A) + diag(stats::runif(n, 0.5, 1.5))
  res  <- waldTest(list(pvals = pv, vcov = vcov),
                   cc       = list(list(coef = treatments,
                                        type = "con",
                                        comp = "pairwise")),
                   test     = "F",
                   df_error = 100L)
  expect_s3_class(plot_waldTest(res), "ggplot")
})

# ---------------------------------------------------------------------------
# 9. p-value labels
# ---------------------------------------------------------------------------
test_that("very small p-value gets 'p<0.001' label", {
  pv   <- data.frame(
    Treatment       = factor(c("A", "B")),
    predicted.value = c(10, 100),
    std.error       = c(0.1, 0.1),
    status          = factor(c("Estimable","Estimable")),
    stringsAsFactors = FALSE
  )
  vcov <- diag(c(0.01, 0.01))
  res  <- waldTest(list(pvals = pv, vcov = vcov),
                   cc = list(list(coef = c("A","B"), type = "con",
                                  comp = c(-1, 1))))
  df   <- plot_waldTest(res, return_data = TRUE)
  expect_equal(df$p_label[1L], "p<0.001")
})

test_that("moderate p-value gets 'p=...' label format", {
  res <- make_wt_result(seed = 200L)
  df  <- plot_waldTest(res, return_data = TRUE)
  mod <- df[df$P.Value >= 0.001, ]
  if (nrow(mod) > 0L)
    expect_true(all(grepl("^p=", mod$p_label)))
})

# ---------------------------------------------------------------------------
# 10. Input validation
# ---------------------------------------------------------------------------
test_that("non-list res errors", {
  expect_error(plot_waldTest("not_a_list"), "'res' must be the list")
})

test_that("res without test/adjust errors", {
  expect_error(plot_waldTest(list(Contrasts = NULL, Zero = NULL)),
               "'res' must be the list")
})

test_that("res with NULL Contrasts and NULL Zero errors", {
  res           <- make_wt_result()
  res$Contrasts <- NULL
  res$Zero      <- NULL
  expect_error(plot_waldTest(res), "no Contrast results")
})

test_that("non-logical facet errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, facet = "yes"), "'facet' must be TRUE or FALSE")
})

test_that("facet of length > 1 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, facet = c(TRUE, FALSE)),
               "'facet' must be TRUE or FALSE")
})

test_that("ci_level = 1 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, ci_level = 1), "'ci_level'")
})

test_that("ci_level = 0 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, ci_level = 0), "'ci_level'")
})

test_that("ci_level > 1 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, ci_level = 1.5), "'ci_level'")
})

test_that("alpha = 0 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, alpha = 0), "'alpha'")
})

test_that("alpha = 1 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, alpha = 1), "'alpha'")
})

test_that("alpha < 0 errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, alpha = -0.1), "'alpha'")
})

test_that("non-theme theme argument errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, theme = "bw"),
               "'theme' must be a ggplot2 theme")
})

test_that("non-logical return_data errors", {
  res <- make_wt_result()
  expect_error(plot_waldTest(res, return_data = "yes"),
               "'return_data' must be TRUE or FALSE")
})

# ---------------------------------------------------------------------------
# 11. Internal helper: .wt_detect_by
# ---------------------------------------------------------------------------
test_that(".wt_detect_by returns NULL when no by column", {
  res <- make_wt_result()
  expect_null(biomAid:::.wt_detect_by(res))
})

test_that(".wt_detect_by detects by column from Contrasts", {
  res <- make_wt_result(by_var = "Site", n_groups = 2L)
  expect_equal(biomAid:::.wt_detect_by(res), "Site")
})

# ---------------------------------------------------------------------------
# 12. Internal helper: .wt_rescale01
# ---------------------------------------------------------------------------
test_that(".wt_rescale01 maps to [0, 1]", {
  x   <- c(2, 4, 6, 8, 10)
  res <- biomAid:::.wt_rescale01(x)
  expect_equal(min(res), 0)
  expect_equal(max(res), 1)
})

test_that(".wt_rescale01 preserves order", {
  x   <- c(1, 5, 3, 9, 2)
  res <- biomAid:::.wt_rescale01(x)
  expect_equal(order(res), order(x))
})

test_that(".wt_rescale01 handles constant vector without error", {
  expect_no_error(biomAid:::.wt_rescale01(rep(5, 4)))
})

test_that(".wt_rescale01 on constant vector returns correct length", {
  res <- biomAid:::.wt_rescale01(rep(3, 3))
  expect_length(res, 3)
})

# ---------------------------------------------------------------------------
# 13. p-value adjustment reflected in data
# ---------------------------------------------------------------------------
test_that("bonferroni adjustment produces p-values >= unadjusted", {
  set.seed(5L)
  treatments <- c("N0","N1","N2")
  n  <- 3L
  pv <- data.frame(
    Treatment       = factor(treatments),
    predicted.value = stats::rnorm(n, 50, 5),
    std.error       = stats::runif(n, 1, 2),
    status          = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(stats::rnorm(n * n) * 0.3, n, n)
  vcov <- crossprod(A) + diag(stats::runif(n, 0.5, 1))
  pred <- list(pvals = pv, vcov = vcov)

  res_none <- waldTest(pred,
                       cc     = list(list(coef = treatments, type = "con",
                                          comp = "pairwise")),
                       adjust = "none")
  res_bonf <- waldTest(pred,
                       cc     = list(list(coef = treatments, type = "con",
                                          comp = "pairwise")),
                       adjust = "bonferroni")

  df_none <- plot_waldTest(res_none, return_data = TRUE)
  df_bonf <- plot_waldTest(res_bonf, return_data = TRUE)
  expect_true(all(df_bonf$P.Value >= df_none$P.Value))
})
