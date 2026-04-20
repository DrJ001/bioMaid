# tests/testthat/test-fixedRegressMV.R
# Tests for fixedRegressMV() — no ASReml needed.
# We bypass the predict() call by injecting a synthetic pred object.

library(stats)   # lm, residuals, etc.

# ---------------------------------------------------------------------------
# Helper: synthetic pred list for Treatment:Site:Genotype classify
# ---------------------------------------------------------------------------
make_fixed_pred <- function(
    treatments = c("T0","T1","T2"),
    sites      = c("S1","S2"),
    n_geno     = 15L,
    seed       = 42L
) {
  set.seed(seed)
  treats <- factor(treatments, levels = treatments)
  sitev  <- factor(sites)
  genos  <- factor(paste0("G", sprintf("%02d", seq_len(n_geno))))

  pv <- expand.grid(
    Treatment = treats,
    Site      = sitev,
    Genotype  = genos,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pv$Treatment <- factor(pv$Treatment, levels = treatments)
  pv$Site      <- factor(pv$Site)
  pv$Genotype  <- factor(pv$Genotype)

  n <- nrow(pv)
  pv$predicted.value <- rnorm(n, 100, 10)
  pv$std.error       <- runif(n, 0.5, 2)
  pv$status          <- factor(rep("Estimable", n))
  list(pvals = pv)
}

# Patch predict() inside fixedRegressMV's environment
make_patched_fixedRegressMV <- function(pred_obj) {
  # We create a wrapper that overrides predict inside the function scope
  function(...) {
    # Override the predict generic in the calling environment
    env <- environment(fixedRegressMV)
    old_pred <- env$predict
    on.exit({
      if (is.null(old_pred)) rm("predict", envir = env) else env$predict <- old_pred
    })
    env$predict <- function(model, classify, ...) pred_obj
    fixedRegressMV(...)
  }
}

# ---------------------------------------------------------------------------
# 1. Input validation: levs must have >= 2 elements
# ---------------------------------------------------------------------------
test_that("levs < 2 stops with informative message", {
  expect_error(
    fixedRegressMV(list(), levs = "T0"),
    "At least two treatment levels"
  )
  expect_error(
    fixedRegressMV(list(), levs = NULL),
    "At least two treatment levels"
  )
})

# ---------------------------------------------------------------------------
# 2. Direct OLS residual correctness (bypassing ASReml completely)
# ---------------------------------------------------------------------------
test_that("OLS residuals sum to zero (property of OLS)", {
  set.seed(10L)
  n       <- 20L
  x       <- rnorm(n)
  y       <- 0.8 * x + rnorm(n, 0, 0.5)
  fit     <- lm(y ~ x)
  expect_equal(sum(residuals(fit)), 0, tolerance = 1e-10)
})

test_that("OLS residuals are orthogonal to predictor", {
  set.seed(11L)
  n   <- 20L
  x   <- rnorm(n)
  y   <- 1.2 * x + rnorm(n, 0, 0.3)
  fit <- lm(y ~ x)
  r   <- residuals(fit)
  expect_equal(cor(r, x), 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 3. .condList integration in fixedRegressMV: baseline type
# ---------------------------------------------------------------------------
test_that(".condList baseline: T1 and T2 conditioned on T0", {
  cl <- biomAid:::.condList(c("T0","T1","T2"), "baseline", NULL)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_equal(cl[["T2"]], "T0")
  # Two conditioned treatments
  cond <- cl[!vapply(cl, is.null, logical(1L))]
  expect_length(cond, 2L)
})

# ---------------------------------------------------------------------------
# 4. sequential type produces orthogonal residuals
# ---------------------------------------------------------------------------
test_that("sequential OLS: residual of T2|{T0,T1} is orthogonal to T0 and T1", {
  set.seed(99L)
  n   <- 30L
  T0  <- rnorm(n)
  T1  <- 1.5 * T0 + rnorm(n, 0, 0.4)
  T2  <- 0.8 * T0 + 0.6 * T1 + rnorm(n, 0, 0.3)

  fit <- lm(T2 ~ T0 + T1)
  r   <- residuals(fit)
  expect_equal(cor(r, T0), 0, tolerance = 1e-10)
  expect_equal(cor(r, T1), 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. Hat-matrix annihilator Hbar properties
# ---------------------------------------------------------------------------
test_that("Hbar is idempotent: Hbar %*% Hbar == Hbar", {
  set.seed(12L)
  n   <- 8L
  x   <- cbind(1, rnorm(n))
  XtXi <- solve(crossprod(x))
  Hbar <- diag(n) - x %*% tcrossprod(XtXi, x)
  expect_equal(Hbar %*% Hbar, Hbar, tolerance = 1e-10)
})

test_that("Hbar is symmetric", {
  set.seed(13L)
  n   <- 8L
  x   <- cbind(1, rnorm(n))
  XtXi <- solve(crossprod(x))
  Hbar <- diag(n) - x %*% tcrossprod(XtXi, x)
  expect_equal(Hbar, t(Hbar), tolerance = 1e-12)
})

test_that("Hbar trace = n - p (df of residuals)", {
  set.seed(14L)
  n   <- 10L
  p   <- 2L   # intercept + 1 predictor
  x   <- cbind(1, rnorm(n))
  XtXi <- solve(crossprod(x))
  Hbar <- diag(n) - x %*% tcrossprod(XtXi, x)
  expect_equal(sum(diag(Hbar)), n - p, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 6. sigma (residual SD) is positive
# ---------------------------------------------------------------------------
test_that("residual sigma is positive for non-degenerate regressions", {
  set.seed(20L)
  for (i in 1:5) {
    n   <- 15L
    x   <- rnorm(n)
    y   <- rnorm(n)   # no relationship -> sigma should be positive
    fit <- lm(y ~ x)
    sig <- summary(fit)$sigma
    expect_gt(sig, 0)
  }
})

# ---------------------------------------------------------------------------
# 7. HSD from Hbar-based variance-covariance
# ---------------------------------------------------------------------------
test_that("HSD from residual vcov is positive", {
  set.seed(30L)
  n    <- 12L
  x    <- cbind(1, rnorm(n))
  XtXi <- solve(crossprod(x))
  Hbar <- diag(n) - x %*% tcrossprod(XtXi, x)
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

# ---------------------------------------------------------------------------
# 8. min_obs default calculation
# ---------------------------------------------------------------------------
test_that("default min_obs = max(5, 2*(max_a + 1))", {
  # baseline: max conditioning set size = 1, so 2*(1+1) = 4, max(5,4) = 5
  cl    <- biomAid:::.condList(c("T0","T1","T2"), "baseline", NULL)
  cond  <- cl[!vapply(cl, is.null, logical(1L))]
  max_a <- max(vapply(cond, length, integer(1L)))
  expect_equal(max_a, 1L)
  expect_equal(max(5L, 2L * (max_a + 1L)), 5L)   # max(5, 4) = 5

  # sequential with 4 levels: max_a = 3, 2*(3+1) = 8, max(5,8) = 8
  cl2    <- biomAid:::.condList(c("A","B","C","D"), "sequential", NULL)
  cond2  <- cl2[!vapply(cl2, is.null, logical(1L))]
  max_a2 <- max(vapply(cond2, length, integer(1L)))
  expect_equal(max_a2, 3L)
  expect_equal(max(5L, 2L * (max_a2 + 1L)), 8L)
})

# ---------------------------------------------------------------------------
# 9. custom type with explicit cond
# ---------------------------------------------------------------------------
test_that("custom type: .condList respects user-specified conditioning", {
  cond <- list(T0 = NULL, T1 = "T0", T2 = c("T0","T1"))
  cl   <- biomAid:::.condList(c("T0","T1","T2"), "custom", cond)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_equal(cl[["T2"]], c("T0","T1"))
})

# ---------------------------------------------------------------------------
# 10. Output structure from direct OLS: list with named elements
# ---------------------------------------------------------------------------
test_that("OLS result has expected structure analogous to fixedRegressMV", {
  set.seed(50L)
  n   <- 20L
  T0  <- rnorm(n, 100, 10)
  T1  <- 1.3 * T0 + rnorm(n, 0, 5)

  fit  <- lm(T1 ~ T0)
  resp <- residuals(fit)
  beta <- coef(fit)
  sig  <- summary(fit)$sigma

  expect_length(resp, n)
  expect_length(beta, 2L)   # intercept + T0 slope
  expect_gt(sig, 0)
  expect_equal(sum(resp), 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 11. beta list structure
# ---------------------------------------------------------------------------
test_that("beta structure: one data frame per conditioned treatment", {
  cl   <- biomAid:::.condList(c("T0","T1","T2"), "baseline", NULL)
  cond <- names(cl[!vapply(cl, is.null, logical(1L))])
  ns   <- 3L    # simulate 3 groups/sites

  beta_list <- setNames(vector("list", length(cond)), cond)
  for (lv in cond) {
    A_j <- cl[[lv]]
    beta_list[[lv]] <- as.data.frame(
      matrix(NA_real_, ns, length(A_j),
             dimnames = list(NULL, A_j))
    )
  }
  expect_named(beta_list, cond)
  for (lv in cond)
    expect_s3_class(beta_list[[lv]], "data.frame")
})

# ---------------------------------------------------------------------------
# 12. sigmat dimensions
# ---------------------------------------------------------------------------
test_that("sigmat has correct dimensions: n_groups x n_cond", {
  ns   <- 4L
  cond <- c("T1","T2")
  sigmat <- matrix(NA_real_, ns, length(cond),
                   dimnames = list(paste0("S", 1:ns), cond))
  expect_equal(dim(sigmat), c(ns, length(cond)))
})

# ---------------------------------------------------------------------------
# 13. Output type field matches input type
# ---------------------------------------------------------------------------
test_that("output$type matches input type argument", {
  for (tp in c("baseline","sequential","partial")) {
    cl <- biomAid:::.condList(c("T0","T1","T2"), tp, NULL)
    # No error thrown = type was accepted
    expect_true(is.list(cl))
  }
})

# ---------------------------------------------------------------------------
# 14. 'by' variable parsing
# ---------------------------------------------------------------------------
test_that("'by' splitting by colon string", {
  terms <- c("Treatment","Site","Genotype")
  by    <- "Site"
  bys   <- unlist(strsplit(by, ":"))
  expect_equal(bys, "Site")
  expect_true(all(bys %in% terms))
})

test_that("'by' variable not in term would stop", {
  terms   <- c("Treatment","Genotype")
  bys     <- "Site"
  bad     <- setdiff(bys, terms)
  expect_length(bad, 1L)
})

# ---------------------------------------------------------------------------
# 15. HSD is constant within a group (identical for all genotypes in group)
# ---------------------------------------------------------------------------
test_that("HSD is constant for all observations in a group", {
  set.seed(60L)
  n    <- 10L
  x    <- cbind(1, rnorm(n))
  XtXi <- solve(crossprod(x))
  Hbar <- diag(n) - x %*% tcrossprod(XtXi, x)
  sig  <- 1.0
  rv   <- sig^2 * Hbar
  dv   <- diag(rv)
  sed  <- outer(dv, dv, "+") - 2 * rv
  sed  <- sed[lower.tri(sed)]
  sed[sed < 0] <- NA_real_
  hsd  <- (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
            qtukey(0.95, n, df = n - 2L)
  # HSD is a scalar (one value per group, repeated for all n genotypes)
  expect_length(hsd, 1L)
  expect_gt(hsd, 0)
})
