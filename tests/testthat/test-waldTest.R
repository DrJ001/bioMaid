# tests/testthat/test-waldTest.R
# Consolidated tests for waldTest() and waldTest.asreml().
# No ASReml licence required: pred objects are built synthetically.
# waldTest.asreml() tests mock predict() via local_mocked_bindings.

# ===========================================================================
# Helper: build a synthetic pred list (pvals + vcov)
# ===========================================================================
make_pred_wt <- function(
    treatments = c("N0","N1","N2"),
    n_per_trt  = 1L,
    seed       = 1L
) {
  set.seed(seed)
  ntreat <- length(treatments)
  if (n_per_trt == 1L) {
    pv <- data.frame(
      Treatment       = factor(treatments, levels = treatments),
      predicted.value = rnorm(ntreat, 50, 5),
      std.error       = runif(ntreat, 0.5, 1.5),
      status          = factor(rep("Estimable", ntreat)),
      stringsAsFactors = FALSE
    )
  } else {
    genotypes <- paste0("G", seq_len(n_per_trt))
    pv <- expand.grid(
      Treatment = factor(treatments, levels = treatments),
      Genotype  = factor(genotypes),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    pv$Treatment       <- factor(pv$Treatment, levels = treatments)
    pv$Genotype        <- factor(pv$Genotype)
    pv$predicted.value <- rnorm(nrow(pv), 50, 5)
    pv$std.error       <- runif(nrow(pv), 0.5, 1.5)
    pv$status          <- factor(rep("Estimable", nrow(pv)))
  }

  n    <- nrow(pv)
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(runif(n, 0.2, 0.8))
  list(pvals = pv, vcov = vcov)
}

# ===========================================================================
# SECTION 1 — Private helpers (accessed via :::)
# ===========================================================================

# ---------------------------------------------------------------------------
# 1a. .rowLabels
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 1b. .pairwiseMat
# ---------------------------------------------------------------------------
test_that(".pairwiseMat: k=2 gives [1,-1]", {
  mat <- biomAid:::.pairwiseMat(2L)
  expect_equal(mat, matrix(c(1,-1), nrow = 1L))
})

test_that(".pairwiseMat: k=3 gives C(3,2)=3 rows, each summing to 0", {
  mat <- biomAid:::.pairwiseMat(3L)
  expect_equal(nrow(mat), 3L)
  expect_equal(ncol(mat), 3L)
  expect_equal(rowSums(mat), rep(0, 3L))
})

test_that(".pairwiseMat: k=4 gives 6 contrasts", {
  mat <- biomAid:::.pairwiseMat(4L)
  expect_equal(nrow(mat), choose(4L, 2L))
})

# ---------------------------------------------------------------------------
# 1c. .resolveCoef
# ---------------------------------------------------------------------------
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

# ===========================================================================
# SECTION 2 — Wald statistic formula
# ===========================================================================

test_that("Wald statistic = estimate^2 / variance", {
  p    <- make_pred_wt(seed = 2L)
  tau  <- p$pvals$predicted.value
  vcov <- p$vcov
  c_i  <- c(-1, 1, 0)
  cv   <- drop(c_i %*% vcov %*% c_i)
  est  <- drop(c_i %*% tau)
  W    <- est^2 / cv
  expect_equal(W, est^2 / cv, tolerance = 1e-12)
})

# ===========================================================================
# SECTION 3 — Input validation
# ===========================================================================

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

test_that("pred with all NA predicted.value errors", {
  p <- make_pred_wt(seed = 13L)
  p$pvals$predicted.value <- NA_real_
  expect_error(
    waldTest(p, cc = list(list(coef = 1L, type = "con", comp = 1L))),
    "No non-missing"
  )
})

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

test_that("by variable not in pred$pvals errors", {
  p <- make_pred_wt(seed = 19L)
  expect_error(
    waldTest(p,
             cc = list(list(coef = c("N0","N1"), type = "con", comp = c(-1,1))),
             by = "NonExistent"),
    "not found in predictions"
  )
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

# ===========================================================================
# SECTION 4 — Contrast ("con") tests
# ===========================================================================

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

test_that("Estimate for (N1-N0) is positive when N1 > N0", {
  p <- make_pred_wt(treatments = c("N0","N1"), seed = 14L)
  p$pvals$predicted.value <- c(10, 20)
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
# NEW TEST 1 — comp as a 1×2 MATRIX (not just vector)
# ---------------------------------------------------------------------------
test_that("comp as a 1x2 matrix works identically to comp as a vector", {
  p    <- make_pred_wt(treatments = c("N0","N1"), seed = 101L)
  vec  <- c(-1, 1)
  mat  <- matrix(vec, nrow = 1L)

  res_vec <- waldTest(p,
                      cc = list(list(coef = c("N0","N1"),
                                     type = "con",
                                     comp = vec)))
  res_mat <- waldTest(p,
                      cc = list(list(coef = c("N0","N1"),
                                     type = "con",
                                     comp = mat)))

  expect_equal(res_mat$Contrasts$Estimate,       res_vec$Contrasts$Estimate)
  expect_equal(res_mat$Contrasts$Wald.Statistic, res_vec$Contrasts$Wald.Statistic)
  expect_equal(res_mat$Contrasts$Comparison,     res_vec$Contrasts$Comparison)
})

# ---------------------------------------------------------------------------
# NEW TEST 2 — comp = NULL triggers a warning (Helmert contrast used)
# ---------------------------------------------------------------------------
test_that("comp = NULL triggers Helmert warning and still returns a result", {
  p <- make_pred_wt(treatments = c("N0","N1","N2"), seed = 102L)
  expect_warning(
    res <- waldTest(p,
                    cc = list(list(coef = c("N0","N1","N2"),
                                   type = "con",
                                   comp = NULL))),
    "Helmert"
  )
  expect_false(is.null(res$Contrasts))
})

# ---------------------------------------------------------------------------
# NEW TEST 3 — group override in .buildCRows changes Comparison label
# ---------------------------------------------------------------------------
test_that("group override changes Comparison label to supplied names", {
  p <- make_pred_wt(treatments = c("N0","N1"), seed = 103L)
  res <- waldTest(p,
                  cc = list(list(coef  = c("N0","N1"),
                                 type  = "con",
                                 comp  = c(-1, 1),
                                 group = list(left  = "Treatment_A",
                                              right = "Treatment_B"))))
  expect_true(grepl("Treatment_A", res$Contrasts$Comparison))
  expect_true(grepl("Treatment_B", res$Contrasts$Comparison))
})

# ---------------------------------------------------------------------------
# NEW TEST 4 — by as a two-element character VECTOR produces pasted column
# ---------------------------------------------------------------------------
test_that("by as two-element vector produces pasted grouping column name", {
  # Build pred with Treatment x Site x Genotype structure so that after
  # excluding Treatment and Site (the `by` variables) there is still a
  # Genotype factor to form within-group row labels.
  set.seed(104L)
  treatments <- c("N0","N1")
  sites      <- c("S1","S2")
  genotypes  <- c("G1","G2")
  pv <- expand.grid(
    Treatment = factor(treatments, levels = treatments),
    Site      = factor(sites),
    Genotype  = factor(genotypes),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pv$Treatment       <- factor(pv$Treatment, levels = treatments)
  pv$Site            <- factor(pv$Site)
  pv$Genotype        <- factor(pv$Genotype)
  pv$predicted.value <- rnorm(nrow(pv), 50, 5)
  pv$std.error       <- runif(nrow(pv), 0.5, 1.5)
  pv$status          <- factor(rep("Estimable", nrow(pv)))
  n <- nrow(pv)
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(runif(n, 0.2, 0.8))
  p <- list(pvals = pv, vcov = vcov)

  # by = c("Site","Treatment") — two-element vector, each group contains
  # 2 Genotype rows; coef refers to within-group Genotype labels
  res <- waldTest(p,
                  cc = list(list(coef = c("G1","G2"),
                                 type = "con",
                                 comp = c(-1, 1))),
                  by = c("Site", "Treatment"))   # two-element by vector

  # The output grouping column should be named "Site:Treatment"
  expect_true("Site:Treatment" %in% names(res$Contrasts))
})

# ===========================================================================
# SECTION 5 — Zero ("zero") equality tests
# ===========================================================================

test_that("waldTest zero test returns Zero data frame", {
  p <- make_pred_wt(seed = 8L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N1","N2"),
                                 type = "zero",
                                 group = "joint")))
  expect_false(is.null(res$Zero))
  expect_equal(res$Zero$df, 2L)
})

test_that("zero test df equals length of coef", {
  p <- make_pred_wt(seed = 20L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1","N2"),
                                 type = "zero",
                                 group = "all_zero")))
  expect_equal(res$Zero$df, 3L)
})

# ===========================================================================
# SECTION 6 — by-group splitting
# ===========================================================================

test_that("by group produces one row per group for single contrast", {
  p <- make_pred_wt(treatments = c("N0","N1"), n_per_trt = 3L, seed = 9L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))),
                  by = "Genotype")
  # 3 genotypes x 1 contrast each = 3 rows
  expect_equal(nrow(res$Contrasts), 3L)
  expect_true("Genotype" %in% names(res$Contrasts))
})

# ===========================================================================
# SECTION 7 — P-value adjustment
# ===========================================================================

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
# NEW TEST 5 — adjust applied when Zero results are present
# ---------------------------------------------------------------------------
test_that("bonferroni adjust is applied to Zero p-values when present", {
  p <- make_pred_wt(treatments = c("N0","N1","N2"), seed = 105L)

  # Without adjustment
  res_none <- waldTest(p,
                       cc = list(
                         list(coef = c("N0"), type = "zero", group = "zero_N0"),
                         list(coef = c("N1"), type = "zero", group = "zero_N1"),
                         list(coef = c("N2"), type = "zero", group = "zero_N2")
                       ),
                       adjust = "none")

  # With Bonferroni adjustment
  res_bonf <- waldTest(p,
                       cc = list(
                         list(coef = c("N0"), type = "zero", group = "zero_N0"),
                         list(coef = c("N1"), type = "zero", group = "zero_N1"),
                         list(coef = c("N2"), type = "zero", group = "zero_N2")
                       ),
                       adjust = "bonferroni")

  raw_pvals  <- res_none$Zero$P.Value
  bonf_pvals <- res_bonf$Zero$P.Value
  # Bonferroni adjusted p-values should be >= raw p-values
  expect_true(all(bonf_pvals >= raw_pvals - 1e-10))
})

# ===========================================================================
# SECTION 8 — Output structure (NULL Contrasts / NULL Zero)
# ===========================================================================

# ---------------------------------------------------------------------------
# NEW TEST 6 — Contrasts-only output: Zero is NULL
# ---------------------------------------------------------------------------
test_that("cc with only 'con' elements gives NULL Zero", {
  p <- make_pred_wt(seed = 106L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_null(res$Zero)
  expect_false(is.null(res$Contrasts))
})

# ---------------------------------------------------------------------------
# NEW TEST 7 — Zero-only output: Contrasts is NULL
# ---------------------------------------------------------------------------
test_that("cc with only 'zero' elements gives NULL Contrasts", {
  p <- make_pred_wt(seed = 107L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "zero",
                                 group = "joint")))
  expect_null(res$Contrasts)
  expect_false(is.null(res$Zero))
})

# ===========================================================================
# SECTION 9 — waldTest.asreml (predict() mocked; no ASReml licence needed)
# ===========================================================================

# ---------------------------------------------------------------------------
# 9a. F-test uses model$nedf
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() F-test uses model$nedf", {
  p <- make_pred_wt(seed = 1L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- structure(list(nedf = 80L), class = "asreml")
  res <- waldTest.asreml(m,
                          classify = "Treatment",
                          cc   = list(list(coef  = c("N0","N1"),
                                           type  = "con",
                                           comp  = c(-1, 1))),
                          test = "F")
  expect_false(is.null(res$Contrasts))
  expect_true("F.Statistic" %in% names(res$Contrasts))
  expect_equal(res$test, "F")
})

# ---------------------------------------------------------------------------
# 9b. Wald test (no df_error needed)
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() Wald test works", {
  p <- make_pred_wt(seed = 2L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- structure(list(nedf = 50L), class = "asreml")
  res <- waldTest.asreml(m,
                          classify = "Treatment",
                          cc   = list(list(coef  = c("N0","N1","N2"),
                                           type  = "con",
                                           comp  = "pairwise")),
                          test = "Wald")
  expect_true("Wald.Statistic" %in% names(res$Contrasts))
  expect_equal(nrow(res$Contrasts), 3L)
})

# ---------------------------------------------------------------------------
# NEW TEST 8 — waldTest.asreml with test="Wald" passes df_error=NULL without error
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() Wald path skips df_error (passes NULL without error)", {
  p <- make_pred_wt(seed = 108L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- structure(list(nedf = 40L), class = "asreml")
  # Wald path: df_error is set to NULL internally; must not error
  expect_no_error(
    waldTest.asreml(m,
                    classify = "Treatment",
                    cc   = list(list(coef  = c("N0","N1"),
                                     type  = "con",
                                     comp  = c(-1, 1))),
                    test = "Wald")
  )
})

# ---------------------------------------------------------------------------
# 9c. zero test
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() zero test works", {
  p <- make_pred_wt(seed = 3L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- structure(list(nedf = 60L), class = "asreml")
  res <- waldTest.asreml(m,
                          classify = "Treatment",
                          cc   = list(list(coef  = c("N1","N2"),
                                           type  = "zero",
                                           group = "joint")),
                          test = "Wald")
  expect_false(is.null(res$Zero))
  expect_equal(res$Zero$df, 2L)
})

# ---------------------------------------------------------------------------
# 9d. errors for non-asreml object
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() errors for non-asreml object", {
  expect_error(
    waldTest.asreml(list(), classify = "Treatment",
                    cc = list(list(coef = "N0", type = "con", comp = 1L))),
    "class.*asreml"
  )
})

# ---------------------------------------------------------------------------
# 9e. errors when nedf is NULL and test = "F"
# ---------------------------------------------------------------------------
test_that("waldTest.asreml() errors when nedf is NULL and test='F'", {
  p <- make_pred_wt(seed = 4L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- structure(list(nedf = NULL), class = "asreml")
  expect_error(
    waldTest.asreml(m, classify = "Treatment",
                    cc   = list(list(coef = c("N0","N1"), type = "con", comp = c(-1,1))),
                    test = "F"),
    "nedf.*NULL"
  )
})

# ---------------------------------------------------------------------------
# 9f. waldTest.asreml is exported and has correct formals
# ---------------------------------------------------------------------------
test_that("waldTest.asreml is a function with expected formals", {
  expect_true(is.function(waldTest.asreml))
  fns <- names(formals(waldTest.asreml))
  expect_true("object"   %in% fns)
  expect_true("classify" %in% fns)
  expect_true("cc"       %in% fns)
})

# ===========================================================================
# SECTION 10 — print.waldTest
# ===========================================================================

# ---------------------------------------------------------------------------
# 10a. Contrasts-only output (no Zero): Contrast header printed
# ---------------------------------------------------------------------------
test_that("print.waldTest with ONLY Contrasts prints contrast header", {
  p <- make_pred_wt(seed = 5L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_output(print(res), "Contrast")
})

# ---------------------------------------------------------------------------
# NEW TEST 9 — Zero-only output: Zero header printed (no Contrasts header)
# ---------------------------------------------------------------------------
test_that("print.waldTest with ONLY Zero results prints Zero header", {
  p <- make_pred_wt(seed = 6L)
  res <- waldTest(p,
                  cc = list(list(coef  = c("N0","N1","N2"),
                                 type  = "zero",
                                 group = "all")))
  out <- paste(capture.output(print(res)), collapse = "\n")
  # "Zero" header should appear
  expect_true(grepl("Zero", out))
  # "Contrast Tests" header should NOT appear (no con elements)
  expect_false(grepl("Contrast Tests", out))
})

# ---------------------------------------------------------------------------
# NEW TEST 10 — Both Contrasts AND Zero: both headers printed
# ---------------------------------------------------------------------------
test_that("print.waldTest with BOTH Contrasts and Zero prints both headers", {
  p <- make_pred_wt(seed = 109L)
  res <- waldTest(p,
                  cc = list(
                    list(coef = c("N0","N1"), type = "con",  comp = c(-1, 1)),
                    list(coef = c("N1","N2"), type = "zero", group = "joint_N1_N2")
                  ))
  out <- capture.output(print(res))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Contrast",    combined))
  expect_true(grepl("Zero",        combined))
})

# ---------------------------------------------------------------------------
# 10b. Adjustment method printed when not "none"
# ---------------------------------------------------------------------------
test_that("print.waldTest() shows p-value adjustment method", {
  p <- make_pred_wt(treatments = c("N0","N1","N2"), seed = 7L)
  res <- waldTest(p,
                  cc     = list(list(coef = c("N0","N1","N2"),
                                     type = "con", comp = "pairwise")),
                  adjust = "bonferroni")
  expect_output(print(res), "bonferroni")
})

# ---------------------------------------------------------------------------
# 10c. print returns invisibly
# ---------------------------------------------------------------------------
test_that("print.waldTest returns invisibly", {
  p <- make_pred_wt(seed = 16L)
  res <- waldTest(p,
                  cc = list(list(coef = c("N0","N1"),
                                 type = "con",
                                 comp = c(-1, 1))))
  expect_invisible(print(res))
})
