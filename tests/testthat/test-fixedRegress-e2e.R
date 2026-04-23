# tests/testthat/test-fixedRegressMV-e2e.R
# End-to-end tests for fixedRegress().
# predict() is imported from stats, so local_mocked_bindings works.
# All OLS calculations (lm, residuals, model.matrix) are standard R
# and run without mocking.

# ---------------------------------------------------------------------------
# Helper: build a synthetic pred object
# ---------------------------------------------------------------------------
make_frm_pred <- function(treatments  = c("T0","T1","T2"),
                          sites       = NULL,
                          n_geno      = 15L,
                          seed        = 1L) {
  set.seed(seed)
  if (is.null(sites)) {
    pv <- expand.grid(
      Treatment       = factor(treatments, levels = treatments),
      Genotype        = factor(paste0("G", sprintf("%02d", seq_len(n_geno)))),
      KEEP.OUT.ATTRS  = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    pv <- expand.grid(
      Treatment       = factor(treatments, levels = treatments),
      Site            = factor(sites),
      Genotype        = factor(paste0("G", sprintf("%02d", seq_len(n_geno)))),
      KEEP.OUT.ATTRS  = FALSE,
      stringsAsFactors = FALSE
    )
    pv$Site <- factor(pv$Site)
  }
  pv$Treatment <- factor(pv$Treatment, levels = treatments)
  pv$Genotype  <- factor(pv$Genotype)
  pv$predicted.value <- rnorm(nrow(pv), 100, 10)
  pv$std.error       <- runif(nrow(pv), 0.5, 2)
  pv$status          <- factor(rep("Estimable", nrow(pv)))
  list(pvals = pv)
}

mock_asreml <- function() structure(list(), class = "asreml")

# ---------------------------------------------------------------------------
# 1. Baseline: returns named list with correct elements
# ---------------------------------------------------------------------------
test_that("fixedRegress() baseline returns correct list structure", {
  p <- make_frm_pred()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"))
  expect_named(res, c("blues","beta","sigmat","cond_list","type"))
  expect_s3_class(res$blues, "data.frame")
  expect_true(is.list(res$beta))
  expect_true(is.matrix(res$sigmat))
  expect_equal(res$type, "baseline")
})

# ---------------------------------------------------------------------------
# 2. Baseline: resp columns and HSD columns present in blues
# ---------------------------------------------------------------------------
test_that("fixedRegress() baseline blues has resp.* and HSD.* columns", {
  p <- make_frm_pred()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"))
  expect_true(all(c("resp.T1","resp.T2") %in% names(res$blues)))
  expect_true(all(c("HSD.T1", "HSD.T2")  %in% names(res$blues)))
  expect_true(all(c("se.T1",  "se.T2")   %in% names(res$blues)))
})

# ---------------------------------------------------------------------------
# 3. Baseline: residuals sum to zero within each group (OLS property)
# ---------------------------------------------------------------------------
test_that("fixedRegress() baseline: resp residuals sum to zero", {
  p <- make_frm_pred(seed = 2L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"))
  expect_equal(sum(res$blues$resp.T1), 0, tolerance = 1e-8)
  expect_equal(sum(res$blues$resp.T2), 0, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# 4. Sequential: residuals orthogonal to T0 (Gram-Schmidt)
# ---------------------------------------------------------------------------
test_that("fixedRegress() sequential: resp.T1 uncorrelated with T0", {
  p <- make_frm_pred(seed = 3L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"),
                         type = "sequential")
  expect_equal(cor(res$blues$resp.T1, res$blues$T0), 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. Sequential: type field matches
# ---------------------------------------------------------------------------
test_that("fixedRegress() sequential: type field is 'sequential'", {
  p <- make_frm_pred(seed = 4L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"),
                         type = "sequential")
  expect_equal(res$type, "sequential")
})

# ---------------------------------------------------------------------------
# 6. Partial: all treatments conditioned
# ---------------------------------------------------------------------------
test_that("fixedRegress() partial: all levs have resp columns", {
  p <- make_frm_pred(seed = 5L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"),
                         type = "partial")
  expect_true(all(c("resp.T0","resp.T1","resp.T2") %in% names(res$blues)))
  expect_equal(res$type, "partial")
})

# ---------------------------------------------------------------------------
# 7. Custom: user-defined conditioning
# ---------------------------------------------------------------------------
test_that("fixedRegress() custom type works", {
  p <- make_frm_pred(seed = 6L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"),
                         type = "custom",
                         cond = list(T0 = NULL, T1 = "T0", T2 = c("T0","T1")))
  expect_equal(res$type, "custom")
  expect_true(all(c("resp.T1","resp.T2") %in% names(res$blues)))
})

# ---------------------------------------------------------------------------
# 8. by argument: separate groups per site
# ---------------------------------------------------------------------------
test_that("fixedRegress() by='Site': separate regressions per site", {
  p <- make_frm_pred(sites = c("S1","S2"), n_geno = 15L, seed = 7L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Site:Genotype",
                         by   = "Site",
                         levs = c("T0","T1","T2"))
  expect_true("Site" %in% names(res$blues))
  expect_true(all(c("S1","S2") %in% res$blues$Site))
})

# ---------------------------------------------------------------------------
# 9. sigmat dimensions: n_groups x n_conditioned
# ---------------------------------------------------------------------------
test_that("fixedRegress() sigmat has correct dimensions", {
  p <- make_frm_pred(sites = c("S1","S2"), n_geno = 15L, seed = 8L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Site:Genotype",
                         by   = "Site",
                         levs = c("T0","T1","T2"))
  expect_equal(dim(res$sigmat), c(2L, 2L))  # 2 sites, 2 conditioned treatments
  expect_true(all(res$sigmat > 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# 10. min_obs: warning when too few common genotypes
# ---------------------------------------------------------------------------
test_that("fixedRegress() warns when common genotypes < min_obs", {
  p <- make_frm_pred(n_geno = 3L, seed = 9L)  # only 3 genotypes < default min_obs=5
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_warning(
    fixedRegress(mock_asreml(),
                    term    = "Treatment:Genotype",
                    levs    = c("T0","T1","T2"),
                    min_obs = 5L),
    "common genotypes"
  )
})

# ---------------------------------------------------------------------------
# 11. HSD values are positive
# ---------------------------------------------------------------------------
test_that("fixedRegress() HSD values are positive", {
  p <- make_frm_pred(seed = 10L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"))
  hsd_vals <- res$blues$HSD.T1
  expect_true(all(!is.na(hsd_vals)))
  expect_true(all(hsd_vals > 0))
})

# ---------------------------------------------------------------------------
# 12. beta list: data frame per conditioned treatment with correct column
# ---------------------------------------------------------------------------
test_that("fixedRegress() beta contains OLS coefficients", {
  p <- make_frm_pred(seed = 11L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  res <- fixedRegress(mock_asreml(),
                         term = "Treatment:Genotype",
                         levs = c("T0","T1","T2"))
  expect_named(res$beta, c("T1","T2"))
  expect_true("T0" %in% names(res$beta$T1))  # T1 conditioned on T0
  expect_true(all(!is.na(res$beta$T1$T0)))
})
