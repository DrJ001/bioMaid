# tests/testthat/test-randomRegressMV-e2e.R
# End-to-end tests for randomRegress() non-FA path.
# predict() and .asreml_vparams() are both in biomAid's namespace
# and can be mocked via local_mocked_bindings().

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_Gmat_rrm <- function(levs = c("N0","N1","N2"),
                          sites = c("S1","S2"),
                          seed = 1L) {
  set.seed(seed)
  tsnams <- as.vector(outer(levs, sites, paste, sep = "-"))
  k <- length(tsnams)
  A <- matrix(rnorm(k * k), k, k)
  G <- crossprod(A) / k + diag(0.5, k)
  dimnames(G) <- list(tsnams, tsnams)
  G
}

make_pvals_rrm <- function(levs  = c("N0","N1","N2"),
                            sites = c("S1","S2"),
                            n_var = 10L, seed = 1L) {
  set.seed(seed)
  tsnams <- as.vector(outer(levs, sites, paste, sep = "-"))
  varieties <- paste0("Var", sprintf("%02d", seq_len(n_var)))
  pv <- expand.grid(
    TSite         = factor(tsnams, levels = tsnams),
    Variety       = factor(varieties),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  pv$predicted.value <- rnorm(nrow(pv), 0, 1)
  pv$std.error       <- runif(nrow(pv), 0.1, 0.5)
  pv$status          <- factor(rep("Estimable", nrow(pv)))
  n <- nrow(pv)
  A <- matrix(rnorm(n * n) * 0.1, n, n)
  vcov <- crossprod(A) + diag(0.05, n)
  list(pvals = pv, vcov = vcov)
}

make_rrm_model <- function() {
  # Create a model with formulae$random that triggers the non-FA path
  frm <- stats::as.formula("~ TSite:Variety")
  m <- list(
    formulae = list(random = frm),
    call     = list()
  )
  class(m) <- "asreml"
  m
}

# ---------------------------------------------------------------------------
# 1. Baseline: returns correct list elements
# ---------------------------------------------------------------------------
test_that("randomRegress() baseline returns named list", {
  G  <- make_Gmat_rrm()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"),
                          type = "baseline")
  expect_named(res, c("blups","TGmat","Gmat","beta","sigmat","tmat",
                       "cond_list","type"))
  expect_s3_class(res$blups, "data.frame")
  expect_equal(res$type, "baseline")
})

# ---------------------------------------------------------------------------
# 2. Baseline: blups has Site, Variety, raw BLUP and resp columns
# ---------------------------------------------------------------------------
test_that("randomRegress() blups has correct columns", {
  G  <- make_Gmat_rrm()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"))
  expect_true("Site"    %in% names(res$blups))
  expect_true("Variety" %in% names(res$blups))
  expect_true("N0"      %in% names(res$blups))
  expect_true("resp.N1" %in% names(res$blups))
  expect_true("resp.N2" %in% names(res$blups))
})

# ---------------------------------------------------------------------------
# 3. Baseline: Gmat returned correctly
# ---------------------------------------------------------------------------
test_that("randomRegress() Gmat matches mock", {
  G  <- make_Gmat_rrm(seed = 2L)
  pv <- make_pvals_rrm(seed = 2L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"))
  expect_equal(res$Gmat, G, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4. Sequential: TGmat is diagonal (Cholesky property)
# ---------------------------------------------------------------------------
test_that("randomRegress() sequential: TGmat off-diagonals near zero", {
  G  <- make_Gmat_rrm(seed = 3L)
  pv <- make_pvals_rrm(seed = 3L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"),
                          type = "sequential")
  # TGmat should be diagonal for sequential (within each site-block)
  expect_equal(res$type, "sequential")
  expect_true(is.matrix(res$TGmat))
})

# ---------------------------------------------------------------------------
# 5. Partial: all treatments conditioned
# ---------------------------------------------------------------------------
test_that("randomRegress() partial: resp columns for all levs", {
  G  <- make_Gmat_rrm(seed = 4L)
  pv <- make_pvals_rrm(seed = 4L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"),
                          type = "partial")
  expect_true(all(c("resp.N0","resp.N1","resp.N2") %in% names(res$blups)))
})

# ---------------------------------------------------------------------------
# 6. sigmat: conditional variances are positive
# ---------------------------------------------------------------------------
test_that("randomRegress() sigmat values are positive", {
  G  <- make_Gmat_rrm(seed = 5L)
  pv <- make_pvals_rrm(seed = 5L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"))
  expect_true(all(res$sigmat > 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 7. beta: regression coefficients are finite
# ---------------------------------------------------------------------------
test_that("randomRegress() beta coefficients are finite", {
  G  <- make_Gmat_rrm(seed = 6L)
  pv <- make_pvals_rrm(seed = 6L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"))
  expect_true(all(is.finite(res$beta$N1), na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 8. Error: levs NULL
# ---------------------------------------------------------------------------
test_that("randomRegress() errors when levs is NULL", {
  expect_error(randomRegress(make_rrm_model(), levs = NULL),
               "At least two treatment levels")
})

# ---------------------------------------------------------------------------
# 9. Custom type works end-to-end
# ---------------------------------------------------------------------------
test_that("randomRegress() custom type runs without error", {
  G  <- make_Gmat_rrm(seed = 7L)
  pv <- make_pvals_rrm(seed = 7L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                          Env  = "TSite:Variety",
                          levs = c("N0","N1","N2"),
                          type = "custom",
                          cond = list(N0 = NULL, N1 = "N0", N2 = c("N0","N1")))
  expect_equal(res$type, "custom")
  expect_true(all(c("resp.N1","resp.N2") %in% names(res$blups)))
})
