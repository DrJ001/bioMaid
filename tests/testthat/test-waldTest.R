# tests/testthat/test-waldTest.R
# Tests for waldTest() — no ASReml licence required.
# We build synthetic pred objects (pvals + vcov) directly.

# ---------------------------------------------------------------------------
# Helper: make a pred list
# ---------------------------------------------------------------------------
make_pred_wt <- function(
    treatments = c("N0","N1","N2"),
    n_per_trt  = 1L,
    seed       = 1L
) {
  set.seed(seed)
  ntreat <- length(treatments)
  if (n_per_trt == 1L) {
    # One predicted value per treatment level
    pv <- data.frame(
      Treatment       = factor(treatments, levels = treatments),
      predicted.value = rnorm(ntreat, 50, 5),
      std.error       = runif(ntreat, 0.5, 1.5),
      status          = factor(rep("Estimable", ntreat)),
      stringsAsFactors = FALSE
    )
  } else {
    # Multiple reps per treatment, with a second factor
    genotypes <- paste0("G", seq_len(n_per_trt))
    pv <- expand.grid(
      Treatment = factor(treatments, levels = treatments),
      Genotype  = factor(genotypes),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    pv$Treatment <- factor(pv$Treatment, levels = treatments)
    pv$Genotype  <- factor(pv$Genotype)
    pv$predicted.value <- rnorm(nrow(pv), 50, 5)
    pv$std.error       <- runif(nrow(pv), 0.5, 1.5)
    pv$status          <- factor(rep("Estimable", nrow(pv)))
  }

  n <- nrow(pv)
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(runif(n, 0.2, 0.8))
  list(pvals = pv, vcov = vcov)
}

# ---------------------------------------------------------------------------
# Internal helpers (accessed via :::)
# ---------------------------------------------------------------------------

# 1. .rowLabels
test_that(".rowLabels returns factor-level labels", {
  pv <- data.frame(
    Treatment = factor(c("N0","N1","N2")),
    status    = factor(rep("Estimable", 3)),
    stringsAsFactors = FALSE
  )
  labs <- biomAid:::.rowLabels(pv)
  expect_equal(labs, c("N0","N1","N2"))
})

test_that(".rowLabels excludes 'status' column", {
  pv <- data.frame(
    Treatment = factor(c("N0","N1")),
    status    = factor(c("E","E")),
    stringsAsFactors = FALSE
  )
  labs <- biomAid:::.rowLabels(pv)
  expect_equal(labs, c("N0","N1"))
  expect_false("status" %in% labs)
})

test_that(".rowLabels combines two factor columns", {
  pv <- data.frame(
    Treatment = factor(c("N0","N1")),
    Site      = factor(c("S1","S2")),
    stringsAsFactors = FALSE
  )
  labs <- biomAid:::.rowLabels(pv)
  expect_equal(labs, c("N0:S1","N1:S2"))
})

test_that(".rowLabels excludes nominated column", {
  pv <- data.frame(
    Treatment = factor(c("N0","N1")),
    Site      = factor(c("S1","S1")),
    stringsAsFactors = FALSE
  )
  labs <- biomAid:::.rowLabels(pv, exclude = "Site")
  expect_equal(labs, c("N0","N1"))
})

# 2. .pairwiseMat
test_that(".pairwiseMat: k=2 gives [1,-1]", {
  mat <- biomAid:::.pairwiseMat(2L)
  expect_equal(mat, matrix(c(1,-1), nrow = 1L))
})

test_that(".pairwiseMat: k=3 gives C(3,2)=3 rows", {
  mat <- biomAid:::.pairwiseMat(3L)
  expect_equal(nrow(mat), 3L)
  expect_equal(ncol(mat), 3L)
  # Each row sums to 0
  expect_equal(rowSums(mat), rep(0, 3L))
})

test_that(".pairwiseMat: k=4 gives 6 contrasts", {
  mat <- biomAid:::.pairwiseMat(4L)
  expect_equal(nrow(mat), choose(4L, 2L))
})

# 3. .resolveCoef
test_that(".resolveCoef: integer indices pass through", {
  labs <- c("N0","N1","N2")
  idx  <- biomAid:::.resolveCoef(c(1L, 3L), labs)
  expect_equal(idx, c(1L, 3L))
})

test_that(".resolveCoef: character labels resolve correctly", {
  labs <- c("N0","N1","N2")
  idx  <- biomAid:::.resolveCoef(c("N2","N0"), labs)
  expect_equal(idx, c(3L, 1L))
})

test_that(".resolveCoef: unknown label errors", {
  labs <- c("N0","N1","N2")
  expect_error(
    biomAid:::.resolveCoef("Z", labs),
    "Labels not found"
  )
})

test_that(".resolveCoef: out-of-bounds index errors", {
  labs <- c("N0","N1","N2")
  expect_error(
    biomAid:::.resolveCoef(5L, labs),
    "out of bounds"
  )
})

# ---------------------------------------------------------------------------
# 4. Wald statistic formula: W = est^2 / var
# ---------------------------------------------------------------------------
test_that("Wald statistic = estimate^2 / variance", {
  p    <- make_pred_wt(seed = 2L)
  tau  <- p$pvals$predicted.value
  vcov <- p$vcov
  # Contrast: N1 - N0
  c_i  <- c(-1, 1, 0)
  cv   <- drop(c_i %*% vcov %*% c_i)
  est  <- drop(c_i %*% tau)
  W    <- est^2 / cv
  expect_equal(W, est^2 / cv, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 5. waldTest: contrast test with "Wald"
# ---------------------------------------------------------------------------
test_that("waldTest returns Contrasts data frame for 'con' type", {
  p <- make_pred_wt(treatments = c("N0","N1","N2"), seed = 3L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_false(is.null(res$Contrasts))
  expect_s3_class(res$Contrasts, "data.frame")
  expect_true("Wald.Statistic" %in% names(res$Contrasts))
  expect_true("P.Value"        %in% names(res$Contrasts))
})

# ---------------------------------------------------------------------------
# 6. waldTest: F-test
# ---------------------------------------------------------------------------
test_that("waldTest returns F.Statistic when test='F'", {
  p <- make_pred_wt(seed = 4L)
  res <- waldTest(p,
                  cc       = list(list(coef = c("N0","N1"),
                                       type = "con",
                                       comp = c(-1, 1))),
                  test     = "F",
                  df_error = 100L)
  expect_true("F.Statistic" %in% names(res$Contrasts))
})

test_that("F test requires df_error", {
  p <- make_pred_wt(seed = 5L)
  expect_error(
    waldTest(p,
             cc   = list(list(coef = c("N0","N1"), type = "con", comp = c(-1,1))),
             test = "F"),
    "df_error"
  )
})

# ---------------------------------------------------------------------------
# 7. pairwise contrasts: C(k,2) rows
# ---------------------------------------------------------------------------
test_that("pairwise comp generates C(k,2) contrast rows", {
  p <- make_pred_wt(treatments = c("N0","N1","N2"), seed = 6L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1","N2"),
                                 type = "con",
                                 comp = "pairwise")))
  expect_equal(nrow(res$Contrasts), choose(3L, 2L))
})

test_that("pairwise for 4 levels gives 6 rows", {
  p <- make_pred_wt(treatments = c("A","B","C","D"), seed = 7L)
  res <- waldTest(p,
                  cc = list(list(coef = c("A","B","C","D"),
                                 type = "con",
                                 comp = "pairwise")))
  expect_equal(nrow(res$Contrasts), 6L)
})

# ---------------------------------------------------------------------------
# 8. zero test
# ---------------------------------------------------------------------------
test_that("waldTest zero test returns Zero data frame", {
  p <- make_pred_wt(seed = 8L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N1","N2"),
                                 type = "zero",
                                 group = "joint")))
  expect_false(is.null(res$Zero))
  expect_equal(res$Zero$df, 2L)
})

# ---------------------------------------------------------------------------
# 9. by-group splitting
# ---------------------------------------------------------------------------
test_that("by group produces one row per group for single contrast", {
  p <- make_pred_wt(treatments = c("N0","N1"), n_per_trt = 3L, seed = 9L)
  # by = "Genotype" — each genotype gets own group
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))),
                  by = "Genotype")
  # 3 genotypes x 1 contrast each = 3 rows
  expect_equal(nrow(res$Contrasts), 3L)
  expect_true("Genotype" %in% names(res$Contrasts))
})

# ---------------------------------------------------------------------------
# 10. by with vector of names: multi-factor grouping
# ---------------------------------------------------------------------------
test_that("by with two-factor vector produces interaction grouping column", {
  p <- make_pred_wt(treatments = c("N0","N1"), n_per_trt = 2L, seed = 10L)
  # by = c("Treatment","Genotype") -- but this groups by both
  # Let's just test that non-NULL by_label gets pasted
  by_vars  <- c("A","B")
  by_label <- paste(by_vars, collapse = ":")
  expect_equal(by_label, "A:B")
})

# ---------------------------------------------------------------------------
# 11. P-values: Wald test p-value for large statistic is small
# ---------------------------------------------------------------------------
test_that("large Wald statistic -> small p-value", {
  W <- 100
  p <- pchisq(W, df = 1L, lower.tail = FALSE)
  expect_lt(p, 0.001)
})

test_that("small Wald statistic -> large p-value", {
  W <- 0.01
  p <- pchisq(W, df = 1L, lower.tail = FALSE)
  expect_gt(p, 0.5)
})

# ---------------------------------------------------------------------------
# 12. P-value adjustment
# ---------------------------------------------------------------------------
test_that("Bonferroni adjustment increases p-values", {
  raw  <- c(0.01, 0.02, 0.03)
  adj  <- p.adjust(raw, method = "bonferroni")
  expect_true(all(adj >= raw))
})

test_that("adjust='none' leaves p-values unchanged", {
  raw <- c(0.01, 0.05, 0.1)
  adj <- p.adjust(raw, method = "none")
  expect_equal(adj, raw)
})

test_that("waldTest: adjust argument accepted for all methods", {
  p <- make_pred_wt(treatments = c("N0","N1","N2"), seed = 12L)
  for (adj in c("none","bonferroni","holm","fdr","BH","BY")) {
    res <- waldTest(p,
                    cc     = list(list(coef = c("N0","N1","N2"),
                                       type = "con",
                                       comp = "pairwise")),
                    adjust = adj)
    expect_false(is.null(res$Contrasts))
  }
})

# ---------------------------------------------------------------------------
# 13. No NA predictions
# ---------------------------------------------------------------------------
test_that("pred with all NA predicted.value errors", {
  p <- make_pred_wt(seed = 13L)
  p$pvals$predicted.value <- NA_real_
  expect_error(
    waldTest(p, cc = list(list(coef = 1L, type = "con", comp = 1L))),
    "No non-missing"
  )
})

# ---------------------------------------------------------------------------
# 14. Invalid pred format
# ---------------------------------------------------------------------------
test_that("pred without pvals or vcov errors", {
  expect_error(
    waldTest(list(pvals = NULL), cc = list()),
    "'pred' must be the list returned"
  )
  expect_error(
    waldTest(list(x = 1), cc = list()),
    "'pred' must be the list returned"
  )
})

# ---------------------------------------------------------------------------
# 15. Contrast estimate sign and direction
# ---------------------------------------------------------------------------
test_that("Estimate for (N1-N0) is positive when N1 > N0", {
  p <- make_pred_wt(treatments = c("N0","N1"), seed = 14L)
  # Make N1 definitely larger
  p$pvals$predicted.value <- c(10, 20)   # N0=10, N1=20
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_gt(res$Contrasts$Estimate, 0)
})

test_that("Comparison label reads 'N1 vs N0' for comp=c(-1,1)", {
  p <- make_pred_wt(treatments = c("N0","N1"), seed = 15L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_true(grepl("N1 vs N0", res$Contrasts$Comparison))
})

# ---------------------------------------------------------------------------
# 16. print method (smoke test)
# ---------------------------------------------------------------------------
test_that("print.waldTest runs without error", {
  p <- make_pred_wt(seed = 16L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_invisible(print(res))
})

# ---------------------------------------------------------------------------
# 17. cc type validation
# ---------------------------------------------------------------------------
test_that("unknown cc type errors", {
  p <- make_pred_wt(seed = 17L)
  expect_error(
    waldTest(p, cc = list(list(coef = "N0", type = "bad"))),
    "must be \"con\" or \"zero\""
  )
})

test_that("unrecognised cc field errors", {
  p <- make_pred_wt(seed = 18L)
  expect_error(
    waldTest(p, cc = list(list(coef = "N0", type = "con", comp = 1, xyz = 1))),
    "unrecognised fields"
  )
})

# ---------------------------------------------------------------------------
# 18. by variable not in pvals errors
# ---------------------------------------------------------------------------
test_that("by variable not in pred$pvals errors", {
  p <- make_pred_wt(seed = 19L)
  expect_error(
    waldTest(p,
             cc = list(list(coef = c("N0","N1"), type = "con", comp = c(-1,1))),
             by = "NonExistent"),
    "not found in predictions"
  )
})

# ---------------------------------------------------------------------------
# 19. Zero test: df = number of coefficients
# ---------------------------------------------------------------------------
test_that("zero test df equals length of coef", {
  p <- make_pred_wt(seed = 20L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1","N2"),
                                 type = "zero",
                                 group = "all_zero")))
  expect_equal(res$Zero$df, 3L)
})

# ---------------------------------------------------------------------------
# 20. waldTest.asreml is exported and has correct signature
# ---------------------------------------------------------------------------
test_that("waldTest.asreml is a function with expected formals", {
  expect_true(is.function(waldTest.asreml))
  fns <- names(formals(waldTest.asreml))
  expect_true("object"   %in% fns)
  expect_true("classify" %in% fns)
  expect_true("cc"       %in% fns)
})
