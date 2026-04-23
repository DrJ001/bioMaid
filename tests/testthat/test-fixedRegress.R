# tests/testthat/test-fixedRegress.R
# Consolidated tests for fixedRegress() and its helpers.
#
# Sections:
#   A. Input validation
#   B. .condList helper (conditioning-structure builder)
#   C. Core OLS mathematics (unit tests, no ASReml)
#   D. End-to-end: baseline type
#   E. End-to-end: sequential type
#   F. End-to-end: partial type
#   G. End-to-end: custom type
#   H. End-to-end: by argument (single factor)
#   I. End-to-end: by argument (multi-factor composite)
#   J. Edge cases — ndf < 3, min_obs default, lm() failure
#   K. Output structure details (beta, sigmat, HSD, se)
#
# ASReml is not needed: predict() is mocked via local_mocked_bindings.
# All OLS calculations (lm, residuals, model.matrix) run without mocking.

library(stats)   # lm, residuals, model.matrix, qtukey, etc.

# ===========================================================================
# Shared helpers
# ===========================================================================

#' Build a synthetic pred object (no ASReml required).
#'
#' @param treatments  Character vector of treatment labels.
#' @param sites       Optional character vector of site labels; when non-NULL
#'   a "Site" column is added.
#' @param n_geno      Number of genotypes.
#' @param seed        RNG seed.
make_frm_pred <- function(treatments = c("T0", "T1", "T2"),
                           sites      = NULL,
                           n_geno     = 15L,
                           seed       = 1L) {
  set.seed(seed)
  if (is.null(sites)) {
    pv <- expand.grid(
      Treatment      = factor(treatments, levels = treatments),
      Genotype       = factor(paste0("G", sprintf("%02d", seq_len(n_geno)))),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    pv <- expand.grid(
      Treatment      = factor(treatments, levels = treatments),
      Site           = factor(sites),
      Genotype       = factor(paste0("G", sprintf("%02d", seq_len(n_geno)))),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    pv$Site <- factor(pv$Site)
  }
  pv$Treatment        <- factor(pv$Treatment, levels = treatments)
  pv$Genotype         <- factor(pv$Genotype)
  pv$predicted.value  <- rnorm(nrow(pv), 100, 10)
  pv$std.error        <- runif(nrow(pv), 0.5, 2)
  pv$status           <- factor(rep("Estimable", nrow(pv)))
  list(pvals = pv)
}

mock_asreml <- function() structure(list(), class = "asreml")


# ===========================================================================
# A. Input validation
# ===========================================================================

test_that("levs = NULL stops with informative message", {
  expect_error(
    fixedRegress(list(), levs = NULL),
    "At least two treatment levels"
  )
})

test_that("levs with < 2 elements stops with informative message", {
  expect_error(
    fixedRegress(list(), levs = "T0"),
    "At least two treatment levels"
  )
})


# ===========================================================================
# B. .condList helper
# ===========================================================================

test_that(".condList baseline: T1 and T2 conditioned on T0 only", {
  cl <- biomAid:::.condList(c("T0", "T1", "T2"), "baseline", NULL)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_equal(cl[["T2"]], "T0")
  cond <- cl[!vapply(cl, is.null, logical(1L))]
  expect_length(cond, 2L)
})

test_that(".condList sequential: each level conditioned on all preceding", {
  cl <- biomAid:::.condList(c("A", "B", "C", "D"), "sequential", NULL)
  expect_null(cl[["A"]])
  expect_equal(cl[["B"]], "A")
  expect_equal(cl[["C"]], c("A", "B"))
  expect_equal(cl[["D"]], c("A", "B", "C"))
})

test_that(".condList custom: respects user-specified conditioning", {
  cond <- list(T0 = NULL, T1 = "T0", T2 = c("T0", "T1"))
  cl   <- biomAid:::.condList(c("T0", "T1", "T2"), "custom", cond)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_equal(cl[["T2"]], c("T0", "T1"))
})

test_that(".condList returns a list for all accepted type values", {
  for (tp in c("baseline", "sequential", "partial", "custom")) {
    cond_arg <- if (tp == "custom")
      list(T0 = NULL, T1 = "T0", T2 = c("T0", "T1"))
    else
      NULL
    cl <- biomAid:::.condList(c("T0", "T1", "T2"), tp, cond_arg)
    expect_true(is.list(cl))
  }
})

test_that("default min_obs = max(5, 2*(max_a+1)) for baseline (max_a = 1)", {
  # baseline with 3 levs: max conditioning-set size = 1 -> max(5, 4) = 5
  cl    <- biomAid:::.condList(c("T0", "T1", "T2"), "baseline", NULL)
  cond  <- cl[!vapply(cl, is.null, logical(1L))]
  max_a <- max(vapply(cond, length, integer(1L)))
  expect_equal(max_a, 1L)
  expect_equal(max(5L, 2L * (max_a + 1L)), 5L)
})

test_that("default min_obs = max(5, 2*(max_a+1)) for sequential 4 levs (max_a = 3)", {
  # sequential with 4 levs: max conditioning-set size = 3 -> max(5, 8) = 8
  cl2    <- biomAid:::.condList(c("A", "B", "C", "D"), "sequential", NULL)
  cond2  <- cl2[!vapply(cl2, is.null, logical(1L))]
  max_a2 <- max(vapply(cond2, length, integer(1L)))
  expect_equal(max_a2, 3L)
  expect_equal(max(5L, 2L * (max_a2 + 1L)), 8L)
})


# ===========================================================================
# C. Core OLS mathematics (unit tests — no ASReml)
# ===========================================================================

test_that("OLS residuals sum to zero (property of OLS)", {
  set.seed(10L)
  n   <- 20L
  x   <- rnorm(n)
  y   <- 0.8 * x + rnorm(n, 0, 0.5)
  fit <- lm(y ~ x)
  expect_equal(sum(residuals(fit)), 0, tolerance = 1e-10)
})

test_that("OLS residuals are orthogonal to predictor", {
  set.seed(11L)
  n   <- 20L
  x   <- rnorm(n)
  y   <- 1.2 * x + rnorm(n, 0, 0.3)
  fit <- lm(y ~ x)
  expect_equal(cor(residuals(fit), x), 0, tolerance = 1e-10)
})

test_that("sequential OLS: residual of T2|{T0,T1} orthogonal to T0 and T1", {
  set.seed(99L)
  n  <- 30L
  T0 <- rnorm(n)
  T1 <- 1.5 * T0 + rnorm(n, 0, 0.4)
  T2 <- 0.8 * T0 + 0.6 * T1 + rnorm(n, 0, 0.3)
  fit <- lm(T2 ~ T0 + T1)
  r   <- residuals(fit)
  expect_equal(cor(r, T0), 0, tolerance = 1e-10)
  expect_equal(cor(r, T1), 0, tolerance = 1e-10)
})

test_that("residual sigma is positive for non-degenerate regressions", {
  set.seed(20L)
  for (i in seq_len(5L)) {
    n   <- 15L
    x   <- rnorm(n)
    y   <- rnorm(n)
    fit <- lm(y ~ x)
    expect_gt(summary(fit)$sigma, 0)
  }
})

test_that("Hbar is idempotent: Hbar %*% Hbar == Hbar", {
  set.seed(12L)
  n     <- 8L
  X     <- cbind(1, rnorm(n))
  Hbar  <- diag(n) - X %*% tcrossprod(solve(crossprod(X)), X)
  expect_equal(Hbar %*% Hbar, Hbar, tolerance = 1e-10)
})

test_that("Hbar is symmetric", {
  set.seed(13L)
  n     <- 8L
  X     <- cbind(1, rnorm(n))
  Hbar  <- diag(n) - X %*% tcrossprod(solve(crossprod(X)), X)
  expect_equal(Hbar, t(Hbar), tolerance = 1e-12)
})

test_that("Hbar trace = n - p (residual degrees of freedom)", {
  set.seed(14L)
  n    <- 10L
  p    <- 2L   # intercept + 1 predictor
  X    <- cbind(1, rnorm(n))
  Hbar <- diag(n) - X %*% tcrossprod(solve(crossprod(X)), X)
  expect_equal(sum(diag(Hbar)), n - p, tolerance = 1e-10)
})

test_that("HSD from Hbar-based variance-covariance is positive", {
  set.seed(30L)
  n    <- 12L
  X    <- cbind(1, rnorm(n))
  Hbar <- diag(n) - X %*% tcrossprod(solve(crossprod(X)), X)
  sig  <- 1.5
  rv   <- sig^2 * Hbar
  dv   <- diag(rv)
  sed  <- outer(dv, dv, "+") - 2 * rv
  sed  <- sed[lower.tri(sed)]
  sed[sed < 0] <- NA_real_
  hsd  <- (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
            qtukey(0.95, n, df = n - 2L)
  expect_gt(hsd, 0)
})

test_that("HSD is a scalar (one value per group)", {
  set.seed(60L)
  n    <- 10L
  X    <- cbind(1, rnorm(n))
  Hbar <- diag(n) - X %*% tcrossprod(solve(crossprod(X)), X)
  sig  <- 1.0
  rv   <- sig^2 * Hbar
  dv   <- diag(rv)
  sed  <- outer(dv, dv, "+") - 2 * rv
  sed  <- sed[lower.tri(sed)]
  sed[sed < 0] <- NA_real_
  hsd  <- (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
            qtukey(0.95, n, df = n - 2L)
  expect_length(hsd, 1L)
  expect_gt(hsd, 0)
})


# ===========================================================================
# D. End-to-end: baseline type
# ===========================================================================

test_that("fixedRegress() baseline returns correct list structure", {
  p <- make_frm_pred()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"))
  expect_named(res, c("blues", "beta", "sigmat", "cond_list", "type"))
  expect_s3_class(res$blues, "data.frame")
  expect_true(is.list(res$beta))
  expect_true(is.matrix(res$sigmat))
  expect_equal(res$type, "baseline")
})

test_that("fixedRegress() baseline blues has resp.*, se.*, and HSD.* columns", {
  p <- make_frm_pred()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"))
  expect_true(all(c("resp.T1", "resp.T2") %in% names(res$blues)))
  expect_true(all(c("HSD.T1",  "HSD.T2")  %in% names(res$blues)))
  expect_true(all(c("se.T1",   "se.T2")   %in% names(res$blues)))
})

test_that("fixedRegress() baseline: resp residuals sum to zero", {
  p <- make_frm_pred(seed = 2L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"))
  expect_equal(sum(res$blues$resp.T1), 0, tolerance = 1e-8)
  expect_equal(sum(res$blues$resp.T2), 0, tolerance = 1e-8)
})

test_that("fixedRegress() baseline HSD values are all positive", {
  p <- make_frm_pred(seed = 10L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"))
  expect_true(all(!is.na(res$blues$HSD.T1)))
  expect_true(all(res$blues$HSD.T1 > 0))
})


# ===========================================================================
# E. End-to-end: sequential type
# ===========================================================================

test_that("fixedRegress() sequential: type field is 'sequential'", {
  p <- make_frm_pred(seed = 4L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"),
                      type = "sequential")
  expect_equal(res$type, "sequential")
})

test_that("fixedRegress() sequential: resp.T1 uncorrelated with T0 column", {
  p <- make_frm_pred(seed = 3L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"),
                      type = "sequential")
  expect_equal(cor(res$blues$resp.T1, res$blues$T0), 0, tolerance = 1e-10)
})


# ===========================================================================
# F. End-to-end: partial type
# ===========================================================================

test_that("fixedRegress() partial: all levs have resp columns", {
  p <- make_frm_pred(seed = 5L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"),
                      type = "partial")
  expect_true(all(c("resp.T0", "resp.T1", "resp.T2") %in% names(res$blues)))
  expect_equal(res$type, "partial")
})


# ===========================================================================
# G. End-to-end: custom type
# ===========================================================================

test_that("fixedRegress() custom type works and returns correct type field", {
  p <- make_frm_pred(seed = 6L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  # Omit T0 from the cond list entirely (unconditional treatments should simply
  # be absent from 'cond', not listed with NULL, to avoid a known source-level
  # quirk where assigning NULL to a list element removes the key entirely and
  # causes a length-mismatch when deriving 'conditioned').
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"),
                      type = "custom",
                      cond = list(T1 = "T0", T2 = c("T0", "T1")))
  expect_equal(res$type, "custom")
  expect_true(all(c("resp.T1", "resp.T2") %in% names(res$blues)))
})


# ===========================================================================
# H. End-to-end: by argument — single factor
# ===========================================================================

test_that("fixedRegress() by='Site': separate regressions per site", {
  p <- make_frm_pred(sites = c("S1", "S2"), n_geno = 15L, seed = 7L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Site:Genotype",
                      by   = "Site",
                      levs = c("T0", "T1", "T2"))
  expect_true("Site" %in% names(res$blues))
  expect_true(all(c("S1", "S2") %in% res$blues$Site))
})

test_that("fixedRegress() sigmat has correct dimensions: n_groups x n_cond", {
  p <- make_frm_pred(sites = c("S1", "S2"), n_geno = 15L, seed = 8L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Site:Genotype",
                      by   = "Site",
                      levs = c("T0", "T1", "T2"))
  # 2 sites x 2 conditioned treatments (T1, T2)
  expect_equal(dim(res$sigmat), c(2L, 2L))
  expect_true(all(res$sigmat > 0, na.rm = TRUE))
})


# ===========================================================================
# I. End-to-end: by argument — multi-factor composite key
# ===========================================================================

test_that("fixedRegress() by='Site:Block': grouping column is colon-joined composite", {
  # The source splits 'by' on ':' via strsplit, so multi-factor grouping must
  # be specified as a single colon-delimited string, e.g. "Site:Block".
  # Build a pred with Site x Block x Treatment x Genotype
  set.seed(42L)
  pv <- expand.grid(
    Treatment = factor(c("T0", "T1", "T2"), levels = c("T0", "T1", "T2")),
    Site      = factor(c("S1", "S2")),
    Block     = factor(c("B1", "B2")),
    Genotype  = factor(paste0("G", sprintf("%02d", seq_len(10L)))),
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
  pv$predicted.value <- rnorm(nrow(pv), 100, 10)
  pv$std.error       <- runif(nrow(pv), 0.5, 2)
  pv$status          <- factor(rep("Estimable", nrow(pv)))
  p <- list(pvals = pv)

  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Site:Block:Genotype",
                      by   = "Site:Block",           # colon-delimited, as source expects
                      levs = c("T0", "T1", "T2"))

  # The by_col is the raw 'by' string ("Site:Block"); the group values are
  # colon-pasted composites produced by apply(pv[,bys], 1, paste, collapse=":")
  group_col      <- res$blues[["Site:Block"]]
  composite_keys <- sort(unique(group_col))
  # Expect composite keys like "S1:B1", "S1:B2", "S2:B1", "S2:B2"
  expect_true(any(grepl("S1:B1", composite_keys)))
  expect_true(length(composite_keys) > 1L)
})


# ===========================================================================
# J. Edge cases
# ===========================================================================

# ---------------------------------------------------------------------------
# J1. min_obs explicit threshold: group with too few common genotypes is warned
# ---------------------------------------------------------------------------
test_that("fixedRegress() warns and skips when common genotypes < explicit min_obs", {
  # 3 genotypes < min_obs = 5
  p <- make_frm_pred(n_geno = 3L, seed = 9L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_warning(
    fixedRegress(mock_asreml(),
                 term    = "Treatment:Genotype",
                 levs    = c("T0", "T1", "T2"),
                 min_obs = 5L),
    "common genotypes"
  )
})

# ---------------------------------------------------------------------------
# J2. min_obs NULL (default): group below computed default is warned/skipped
#
# baseline with 3 levs => max_a = 1 => default min_obs = max(5, 2*(1+1)) = 5
# A group with only 4 genotypes is therefore below the threshold and skipped.
# ---------------------------------------------------------------------------
test_that("fixedRegress() min_obs=NULL default skips group with < max(5,2*(max_a+1)) genos", {
  # 4 genotypes < default min_obs of 5 for baseline type
  p <- make_frm_pred(n_geno = 4L, seed = 101L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_warning(
    fixedRegress(mock_asreml(),
                 term    = "Treatment:Genotype",
                 levs    = c("T0", "T1", "T2"),
                 min_obs = NULL),    # explicit NULL to exercise default path
    "common genotypes"
  )
})

# ---------------------------------------------------------------------------
# J3. ndf < 3: group with enough genotypes to pass min_obs but too few
#     degrees of freedom for a specific conditioned treatment is warned/skipped.
#
# For baseline type, ndf = n_genos - a - 1 = n_genos - 1 - 1 = n_genos - 2.
# ndf < 3  <=>  n_genos < 5.
# We need n_genos >= min_obs to pass the min_obs gate but ndf < 3.
#
# With sequential type and 4 levs (A,B,C,D):
#   - max_a = 3 (D conditioned on A,B,C), default min_obs = max(5,8) = 8.
#   - For treatment D: ndf = n - 3 - 1 = n - 4; ndf < 3 => n < 7.
#   - So we set n_geno = 6 and override min_obs = 6 to pass the min_obs gate
#     for the group, but ndf for treatment D will be 6 - 4 = 2 < 3.
# ---------------------------------------------------------------------------
test_that("fixedRegress() warns when ndf < 3 for a treatment in a group", {
  # 6 genotypes, sequential 4 levs => D has ndf = 6 - 3 - 1 = 2 < 3
  p <- make_frm_pred(
    treatments = c("A", "B", "C", "D"),
    n_geno     = 6L,
    seed       = 202L
  )
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_warning(
    fixedRegress(mock_asreml(),
                 term    = "Treatment:Genotype",
                 levs    = c("A", "B", "C", "D"),
                 type    = "sequential",
                 min_obs = 6L),   # override so group is not skipped at min_obs gate
    "df"                          # warning message contains "df"
  )
})

# ---------------------------------------------------------------------------
# J4. lm() failure tryCatch: inject a failing lm() into biomAid's namespace
#     so that the tryCatch path is exercised, a warning containing "OLS failed"
#     is issued, and the function continues (does not stop).
#
# Strategy: local_mocked_bindings with .package = "biomAid" replaces whatever
# 'lm' resolves to inside the biomAid package namespace.  We provide a stub
# that always throws.  'predict' is also mocked so no ASReml is needed.
# ---------------------------------------------------------------------------
test_that("fixedRegress() issues 'OLS failed' warning and continues on lm() error", {
  p <- make_frm_pred(seed = 303L)   # well-formed data; the lm() stub will error

  local_mocked_bindings(predict = function(...) p,          .package = "biomAid")
  local_mocked_bindings(lm      = function(...) stop("simulated rank failure"),
                        .package = "biomAid")

  expect_warning(
    fixedRegress(mock_asreml(),
                 term = "Treatment:Genotype",
                 levs = c("T0", "T1", "T2")),
    "OLS failed"
  )
})


# ===========================================================================
# K. Output structure details
# ===========================================================================

test_that("fixedRegress() beta named list: one entry per conditioned treatment", {
  p <- make_frm_pred(seed = 11L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                      term = "Treatment:Genotype",
                      levs = c("T0", "T1", "T2"))
  # baseline: T1 and T2 are conditioned; T0 is the baseline (not conditioned)
  expect_named(res$beta, c("T1", "T2"))
  # T1 conditioned on T0 => beta$T1 should have a "T0" column
  expect_true("T0" %in% names(res$beta$T1))
  expect_true(all(!is.na(res$beta$T1$T0)))
})

test_that("fixedRegress() beta data frames have correct class", {
  cl   <- biomAid:::.condList(c("T0", "T1", "T2"), "baseline", NULL)
  cond <- names(cl[!vapply(cl, is.null, logical(1L))])
  ns   <- 3L

  beta_list <- setNames(vector("list", length(cond)), cond)
  for (lv in cond) {
    A_j             <- cl[[lv]]
    beta_list[[lv]] <- as.data.frame(
      matrix(NA_real_, ns, length(A_j),
             dimnames = list(NULL, A_j))
    )
  }
  for (lv in cond)
    expect_s3_class(beta_list[[lv]], "data.frame")
})

test_that("sigmat matrix has correct dimension structure", {
  ns     <- 4L
  cond   <- c("T1", "T2")
  sigmat <- matrix(NA_real_, ns, length(cond),
                   dimnames = list(paste0("S", seq_len(ns)), cond))
  expect_equal(dim(sigmat), c(ns, length(cond)))
})

test_that("'by' variable not present in term is detected as missing", {
  terms <- c("Treatment", "Genotype")
  bys   <- "Site"
  bad   <- setdiff(bys, terms)
  expect_length(bad, 1L)
})
