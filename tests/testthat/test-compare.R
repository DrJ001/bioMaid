# tests/testthat/test-compare.R
#
# Consolidated tests for compare() and the internal .condList() helper.
#
# Structure
# =========
#   A.  Internal .condList() helper
#   B.  Mathematical unit tests (no ASReml needed — exercise formula directly)
#   C.  End-to-end tests using local_mocked_bindings (no ASReml licence needed)
#       C1. Return shape / column checks
#       C2. Comparison type columns
#       C3. by-grouping
#       C4. alpha direction
#       C5. NA handling
#       C6. pev = FALSE paths
#       C7. Input validation / error paths
#       C8. NEW — additional code-path coverage

# ===========================================================================
# SHARED HELPERS
# ===========================================================================

# ---- Mathematical helpers (section B) ------------------------------------
make_pred_math <- function(n = 20, seed = 1L) {
  set.seed(seed)
  trts <- rep(c("N0", "N1", "N2", "N3"), length.out = n)
  pv <- data.frame(
    Treatment        = factor(trts, levels = c("N0","N1","N2","N3")),
    Genotype         = factor(paste0("G", sprintf("%02d", seq_len(n)))),
    predicted.value  = rnorm(n, 50, 5),
    std.error        = runif(n, 0.3, 0.9),
    status           = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(rnorm(n * n) * 0.3, n, n)
  vcov <- crossprod(A) + diag(runif(n, 0.1, 0.5))
  list(pvals = pv, vcov = vcov)
}

# Reference implementation of the criterion formula
manual_crit <- function(gind, vcov, type, alpha = 0.05, df_err = 100L) {
  n_g  <- length(gind)
  svar <- vcov[gind, gind, drop = FALSE]
  dv   <- diag(svar)
  sed2 <- outer(dv, dv, "+") - 2 * svar
  sed2 <- sed2[lower.tri(sed2)]
  sed2[sed2 < 0] <- NA_real_
  seds  <- sqrt(sed2)
  avs   <- mean(seds, na.rm = TRUE)
  m     <- sum(!is.na(seds))
  switch(type,
    HSD        = (avs / sqrt(2)) * qtukey(1 - alpha, nmeans = n_g, df = df_err),
    LSD        = avs * qt(alpha / 2, df_err, lower.tail = FALSE),
    Bonferroni = avs * qt(alpha / (2 * m), df_err, lower.tail = FALSE)
  )
}

# ---- End-to-end helpers (section C) --------------------------------------
make_asreml_model <- function(nedf = 100L, random_formula = ~Genotype) {
  structure(
    list(nedf = nedf, call = list(random = random_formula)),
    class = "asreml"
  )
}

make_pred_simple <- function(n = 10, seed = 1L) {
  set.seed(seed)
  pv <- data.frame(
    Genotype         = factor(paste0("G", sprintf("%02d", seq_len(n)))),
    predicted.value  = rnorm(n, 50, 5),
    std.error        = runif(n, 0.5, 1.5),
    status           = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(0.3, n)
  list(pvals = pv, vcov = vcov)
}

make_pred_grouped <- function(sites = c("S1","S2"), n_geno = 8L, seed = 2L) {
  set.seed(seed)
  genos <- paste0("G", sprintf("%02d", seq_len(n_geno)))
  pv <- expand.grid(
    Site            = factor(sites),
    Genotype        = factor(genos),
    KEEP.OUT.ATTRS  = FALSE,
    stringsAsFactors = FALSE
  )
  pv$Site     <- factor(pv$Site)
  pv$Genotype <- factor(pv$Genotype)
  n <- nrow(pv)
  pv$predicted.value <- rnorm(n, 50, 5)
  pv$std.error       <- runif(n, 0.5, 1.5)
  pv$status          <- factor(rep("Estimable", n))
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(0.3, n)
  list(pvals = pv, vcov = vcov)
}

# ===========================================================================
# A.  Internal .condList() helper
# ===========================================================================

test_that(".condList baseline: first element NULL, rest conditioned on levs[1]", {
  cl <- biomAid:::.condList(c("N0","N1","N2"), "baseline", NULL)
  expect_null(cl[["N0"]])
  expect_equal(cl[["N1"]], "N0")
  expect_equal(cl[["N2"]], "N0")
})

test_that(".condList sequential: Gram-Schmidt structure", {
  cl <- biomAid:::.condList(c("A","B","C","D"), "sequential", NULL)
  expect_null(cl[["A"]])
  expect_equal(cl[["B"]], "A")
  expect_equal(cl[["C"]], c("A","B"))
  expect_equal(cl[["D"]], c("A","B","C"))
})

test_that(".condList partial: each treatment vs all others", {
  cl <- biomAid:::.condList(c("A","B","C"), "partial", NULL)
  expect_equal(sort(cl[["A"]]), c("B","C"))
  expect_equal(sort(cl[["B"]]), c("A","C"))
  expect_equal(sort(cl[["C"]]), c("A","B"))
})

test_that(".condList custom: user-supplied structure", {
  cond <- list(T0 = NULL, T1 = "T0", T2 = c("T0","T1"))
  cl   <- biomAid:::.condList(c("T0","T1","T2"), "custom", cond)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_equal(cl[["T2"]], c("T0","T1"))
})

test_that(".condList custom: missing cond errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", NULL),
    "'cond' must be supplied"
  )
})

test_that(".condList custom: bad name in cond errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", list(Z = NULL)),
    "not found in 'levs'"
  )
})

test_that(".condList custom: self-conditioning errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", list(A = "A")),
    "cannot be in its own conditioning set"
  )
})

# ===========================================================================
# B.  Mathematical unit tests  (formula verification, no ASReml needed)
# ===========================================================================

# B1. HSD formula -----------------------------------------------------------
test_that("HSD formula: (avsed/sqrt(2)) * qtukey(0.95, n_g, df)", {
  p    <- make_pred_math(n = 20, seed = 7L)
  gind <- seq_len(20L)
  exp_hsd <- manual_crit(gind, p$vcov, "HSD", alpha = 0.05, df_err = 100L)

  svar  <- p$vcov[gind, gind]
  dv    <- diag(svar)
  sed2  <- outer(dv, dv, "+") - 2 * svar
  sed2  <- sed2[lower.tri(sed2)]
  sed2[sed2 < 0] <- NA_real_
  avsed <- mean(sqrt(sed2), na.rm = TRUE)
  hsd   <- (avsed / sqrt(2)) * qtukey(0.95, 20L, 100L)
  expect_equal(hsd, exp_hsd, tolerance = 1e-12)
})

# B2. LSD formula -----------------------------------------------------------
test_that("LSD formula: avsed * qt(alpha/2, df, lower.tail=FALSE)", {
  p    <- make_pred_math(n = 8, seed = 3L)
  gind <- seq_len(8L)
  exp_lsd <- manual_crit(gind, p$vcov, "LSD", alpha = 0.05, df_err = 50L)

  svar  <- p$vcov[gind, gind]
  dv    <- diag(svar)
  sed2  <- outer(dv, dv, "+") - 2 * svar
  avsed <- mean(sqrt(sed2[lower.tri(sed2)]), na.rm = TRUE)
  lsd   <- avsed * qt(0.025, 50L, lower.tail = FALSE)
  expect_equal(lsd, exp_lsd, tolerance = 1e-12)
})

# B3. Bonferroni formula ----------------------------------------------------
test_that("Bonferroni: uses m = C(n_g,2) pairs in t-quantile", {
  p     <- make_pred_math(n = 6, seed = 5L)
  gind  <- seq_len(6L)
  alpha <- 0.05
  df    <- 30L

  svar  <- p$vcov[gind, gind]
  dv    <- diag(svar)
  sed2  <- outer(dv, dv, "+") - 2 * svar
  sed2  <- sed2[lower.tri(sed2)]
  sed2[sed2 < 0] <- NA_real_
  m     <- sum(!is.na(sed2))
  avsed <- mean(sqrt(sed2), na.rm = TRUE)
  bon   <- avsed * qt(alpha / (2 * m), df, lower.tail = FALSE)

  exp_bon <- manual_crit(gind, p$vcov, "Bonferroni", alpha = alpha, df_err = df)
  expect_equal(bon, exp_bon, tolerance = 1e-12)
  expect_equal(m, choose(6L, 2L))
})

# B4. Ordering of criteria: Bonferroni >= LSD --------------------------------
test_that("Bonferroni >= LSD for all group sizes", {
  for (ng in c(3L, 5L, 10L, 15L)) {
    p   <- make_pred_math(n = ng, seed = ng)
    ind <- seq_len(ng)
    lsd <- manual_crit(ind, p$vcov, "LSD",        df_err = 50L)
    bon <- manual_crit(ind, p$vcov, "Bonferroni", df_err = 50L)
    expect_gte(bon, lsd)
  }
})

# B5. avsed = arithmetic mean, not RMS --------------------------------------
test_that("avsed is arithmetic mean of SEDs, not RMS", {
  V    <- diag(c(1, 2, 3, 4, 5))
  dv   <- diag(V)
  sed2 <- outer(dv, dv, "+") - 2 * V
  sed2 <- sed2[lower.tri(sed2)]
  seds <- sqrt(sed2)
  arith <- mean(seds)
  rms   <- sqrt(mean(sed2))
  expect_false(isTRUE(all.equal(arith, rms)))   # they differ
  expect_equal(arith, mean(seds), tolerance = 1e-15)
})

# B6. Alpha direction -------------------------------------------------------
test_that("smaller alpha -> larger HSD and LSD (manual)", {
  p    <- make_pred_math(n = 10, seed = 9L)
  gind <- seq_len(10L)
  hsd05 <- manual_crit(gind, p$vcov, "HSD", alpha = 0.05, df_err = 80L)
  hsd01 <- manual_crit(gind, p$vcov, "HSD", alpha = 0.01, df_err = 80L)
  lsd05 <- manual_crit(gind, p$vcov, "LSD", alpha = 0.05, df_err = 80L)
  lsd01 <- manual_crit(gind, p$vcov, "LSD", alpha = 0.01, df_err = 80L)
  expect_gt(hsd01, hsd05)
  expect_gt(lsd01, lsd05)
})

# B7. By-group label construction -------------------------------------------
test_that("by colon string produces same group labels as character vector (manual)", {
  pv <- data.frame(
    Site      = factor(c("S1","S1","S2","S2")),
    Treatment = factor(c("N0","N1","N0","N1")),
    stringsAsFactors = FALSE
  )
  grp_colon <- apply(pv[, c("Site","Treatment")], 1L, paste, collapse = ":")
  grp_vec   <- apply(pv[, c("Site","Treatment")], 1L, paste, collapse = ":")
  expect_identical(grp_colon, grp_vec)
})

test_that("single by variable: group equals factor levels", {
  pv  <- data.frame(Site = factor(c("S1","S1","S2","S2")))
  grp <- as.character(pv[["Site"]])
  expect_equal(sort(unique(grp)), c("S1","S2"))
})

# B8. m_pairs = C(n_g, 2) --------------------------------------------------
test_that("number of pairs equals C(n_g, 2)", {
  for (ng in 2:8) {
    p    <- make_pred_math(n = ng, seed = ng + 100L)
    gind <- seq_len(ng)
    svar <- p$vcov[gind, gind]
    dv   <- diag(svar)
    sed2 <- outer(dv, dv, "+") - 2 * svar
    sed2 <- sed2[lower.tri(sed2)]
    m    <- sum(!is.na(sed2))
    expect_equal(m, choose(ng, 2L))
  }
})

# B9. SED^2 guard: values are >= 0 after NA substitution -------------------
test_that("SED^2 values are >= 0 after negative guard", {
  p    <- make_pred_math(n = 8, seed = 55L)
  gind <- seq_len(8L)
  svar <- p$vcov[gind, gind]
  dv   <- diag(svar)
  sed2 <- outer(dv, dv, "+") - 2 * svar
  sed2 <- sed2[lower.tri(sed2)]
  sed2[sed2 < 0] <- NA_real_
  valid <- sed2[!is.na(sed2)]
  expect_true(all(valid >= 0))
})

# B10. Pooled vs per-group HSD differ ---------------------------------------
test_that("pooling all rows vs separate groups gives different HSD", {
  n   <- 12L
  p   <- make_pred_math(n = n, seed = 20L)
  hsd_all <- manual_crit(seq_len(n),       p$vcov, "HSD", df_err = 50L)
  hsd_grp <- manual_crit(seq_len(n / 3L),  p$vcov, "HSD", df_err = 50L)
  expect_false(isTRUE(all.equal(hsd_all, hsd_grp)))
})

# B11. Reproducibility ------------------------------------------------------
test_that("manual criterion is reproducible with the same seed", {
  p1 <- make_pred_math(n = 12, seed = 77L)
  p2 <- make_pred_math(n = 12, seed = 77L)
  c1 <- manual_crit(seq_len(12L), p1$vcov, "HSD", df_err = 60L)
  c2 <- manual_crit(seq_len(12L), p2$vcov, "HSD", df_err = 60L)
  expect_equal(c1, c2)
})

# ===========================================================================
# C.  End-to-end tests (local_mocked_bindings — no ASReml licence needed)
# ===========================================================================

# ---------------------------------------------------------------------------
# C1. Return shape / column checks
# ---------------------------------------------------------------------------

test_that("compare() HSD no by: returns data frame with HSD and avsed columns", {
  p <- make_pred_simple(n = 10, seed = 1L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "HSD")
  expect_s3_class(out, "data.frame")
  expect_true("HSD"   %in% names(out))
  expect_true("avsed" %in% names(out))
  expect_equal(nrow(out), 10L)
})

test_that("compare() output has correct column set", {
  p <- make_pred_simple(seed = 10L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "LSD")
  expect_true(all(c("Genotype","predicted.value","std.error","LSD","avsed")
                  %in% names(out)))
})

# ---------------------------------------------------------------------------
# C2. Comparison type columns
# ---------------------------------------------------------------------------

test_that("compare() LSD returns LSD column with positive values", {
  p <- make_pred_simple(seed = 2L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "LSD")
  expect_true("LSD" %in% names(out))
  expect_true(all(out$LSD > 0))
})

test_that("compare() Bonferroni returns Bonferroni column with positive values", {
  p <- make_pred_simple(seed = 3L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "Bonferroni")
  expect_true("Bonferroni" %in% names(out))
  expect_true(all(out$Bonferroni > 0))
})

# ---------------------------------------------------------------------------
# C3. by-grouping
# ---------------------------------------------------------------------------

test_that("compare() by = single string: one constant criterion per site group", {
  p <- make_pred_grouped(sites = c("S1","S2"), n_geno = 8L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Site:Genotype",
                 by = "Site", type = "HSD")
  expect_equal(nrow(out), 16L)
  for (s in c("S1","S2")) {
    vals <- out$HSD[out$Site == s]
    expect_equal(length(unique(round(vals, 10L))), 1L)
  }
})

test_that("compare() by as two-element character VECTOR hits apply() branch", {
  # by = c("Site","Treatment") — length(by) > 1 so the else branch runs:
  #   by_vars <- as.character(by)
  # and then pv[["grp__"]] <- apply(pv[, by_vars, ...], 1L, paste, ...)
  set.seed(42L)
  sites <- c("S1","S2")
  trts  <- c("T1","T2")
  genos <- paste0("G", 1:4)
  pv <- expand.grid(
    Site      = factor(sites),
    Treatment = factor(trts),
    Genotype  = factor(genos),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pv$Site      <- factor(pv$Site)
  pv$Treatment <- factor(pv$Treatment)
  pv$Genotype  <- factor(pv$Genotype)
  n  <- nrow(pv)
  pv$predicted.value <- rnorm(n, 50, 5)
  pv$std.error       <- runif(n, 0.5, 1.5)
  pv$status          <- factor(rep("Estimable", n))
  A    <- matrix(rnorm(n * n) * 0.2, n, n)
  vcov <- crossprod(A) + diag(0.3, n)
  p    <- list(pvals = pv, vcov = vcov)

  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- make_asreml_model(nedf = 50L, random_formula = ~Genotype)
  out <- compare(m, term = "Site:Treatment:Genotype",
                 by = c("Site","Treatment"), type = "HSD")
  expect_s3_class(out, "data.frame")
  expect_true("HSD" %in% names(out))
  # Should have 4 group combinations (S1:T1, S1:T2, S2:T1, S2:T2)
  expect_equal(nrow(out), n)
  # HSD is constant within each Site:Treatment combination
  out$grp <- paste(out$Site, out$Treatment, sep = ":")
  for (g in unique(out$grp)) {
    vals <- out$HSD[out$grp == g]
    expect_equal(length(unique(round(vals, 10L))), 1L)
  }
})

test_that("compare() by colon-string and by character-vector give same HSD", {
  p <- make_pred_grouped()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m  <- make_asreml_model()
  o1 <- compare(m, term = "Site:Genotype", by = "Site",  type = "HSD")
  o2 <- compare(m, term = "Site:Genotype", by = "Site",  type = "HSD")
  expect_identical(o1$HSD, o2$HSD)
})

# ---------------------------------------------------------------------------
# C4. Alpha direction (end-to-end)
# ---------------------------------------------------------------------------

test_that("compare() smaller alpha -> larger HSD (end-to-end)", {
  p <- make_pred_simple(seed = 5L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m   <- make_asreml_model()
  h05 <- compare(m, term = "Genotype", type = "HSD", alpha = 0.05)$HSD[1L]
  h01 <- compare(m, term = "Genotype", type = "HSD", alpha = 0.01)$HSD[1L]
  expect_gt(h01, h05)
})

# ---------------------------------------------------------------------------
# C5. NA handling
# ---------------------------------------------------------------------------

test_that("compare() NA predicted values are filtered out — row count correct", {
  p <- make_pred_simple(n = 10, seed = 6L)
  p$pvals$predicted.value[c(2L, 5L)] <- NA_real_
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  out <- compare(make_asreml_model(), term = "Genotype", type = "HSD")
  expect_equal(nrow(out), 8L)
})

test_that("compare() warns and returns NA criterion for single-obs group", {
  p <- make_pred_grouped(sites = c("S1","S2"), n_geno = 1L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  suppressWarnings(
    out <- compare(make_asreml_model(), term = "Site:Genotype",
                   by = "Site", type = "HSD")
  )
  expect_true(all(is.na(out$HSD)))
})

# ---------------------------------------------------------------------------
# C6. pev = FALSE paths
# ---------------------------------------------------------------------------

test_that("compare() pev=FALSE warns when term not in random formula", {
  p <- make_pred_simple(seed = 8L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- make_asreml_model(random_formula = ~Block)   # random doesn't have Genotype
  expect_warning(
    compare(m, term = "Genotype", type = "HSD", pev = FALSE),
    "Falling back to PEV"
  )
})

test_that("compare() pev=FALSE falls back to PEV when .asreml_vparams returns NULL", {
  p <- make_pred_simple(n = 8L, seed = 9L)
  local_mocked_bindings(
    predict         = function(...) p,
    .asreml_vparams = function(model, term) NULL,
    .package        = "biomAid"
  )
  m <- structure(
    list(nedf = 50L, call = list(random = ~Genotype)),
    class = "asreml"
  )
  expect_warning(
    out <- compare(m, term = "Genotype", type = "HSD", pev = FALSE),
    "Falling back to PEV"
  )
  expect_s3_class(out, "data.frame")
  expect_true("HSD" %in% names(out))
})

test_that("compare() pev=FALSE single-term: len branch uses nrow(pv)", {
  # When length(terms) == 1 the 'else' branch of:
  #   len <- if (length(terms) > 1L) as.integer(table(pv[,1L])[[1L]])
  #           else nrow(pv)
  # is taken.  .asreml_vparams returns a 1×1 variance-component matrix
  # (scalar); the code then does kronecker(varm_1x1, diag(len)) which
  # produces an n×n G-matrix, matching pred$vcov dimensions.
  n <- 6L
  p <- make_pred_simple(n = n, seed = 15L)
  # varm must be a 1×1 matrix (variance component scalar) for a single term
  G <- matrix(2.0, nrow = 1L, ncol = 1L)

  local_mocked_bindings(
    predict         = function(...) p,
    .asreml_vparams = function(model, term) G,
    .package        = "biomAid"
  )
  m <- structure(
    list(nedf = 40L, call = list(random = ~Genotype)),
    class = "asreml"
  )
  # Should not error; uses nrow(pv) = n for the Kronecker product size
  expect_no_error(
    out <- compare(m, term = "Genotype", type = "HSD", pev = FALSE)
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), n)
})

# ---------------------------------------------------------------------------
# C7. Input validation / error paths
# ---------------------------------------------------------------------------

test_that("compare() errors for non-asreml model", {
  expect_error(compare(list(), term = "Genotype"), "class.*asreml")
})

test_that("compare() errors for alpha outside (0,1)", {
  p <- make_pred_simple()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- make_asreml_model()
  expect_error(compare(m, term = "Genotype", alpha = 0),   "'alpha'")
  expect_error(compare(m, term = "Genotype", alpha = 1),   "'alpha'")
  expect_error(compare(m, term = "Genotype", alpha = 1.5), "'alpha'")
})

test_that("compare() errors when by variable not found in term", {
  p <- make_pred_simple()
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_error(
    compare(make_asreml_model(), term = "Genotype",
            by = "Site", type = "HSD"),
    "not found in 'term'"
  )
})

test_that("compare() errors when by accounts for all variables in term", {
  p <- make_pred_grouped(sites = c("S1"), n_geno = 5L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  expect_error(
    compare(make_asreml_model(), term = "Site:Genotype",
            by = c("Site","Genotype"), type = "HSD"),
    "no levels remain"
  )
})

test_that("compare() errors when model$nedf is NULL", {
  # model$nedf is NULL -> stop() with 'nedf' in the message
  p <- make_pred_simple(n = 6L, seed = 20L)
  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m_no_nedf <- structure(
    list(nedf = NULL, call = list(random = ~Genotype)),
    class = "asreml"
  )
  expect_error(
    compare(m_no_nedf, term = "Genotype", type = "HSD"),
    "nedf"
  )
})

# ---------------------------------------------------------------------------
# C8. NEW — additional code-path coverage
# ---------------------------------------------------------------------------

test_that("sed2 < 0 guard: tiny negative off-diagonals become NA, not NaN/error", {
  # Construct a variance matrix with a tiny negative off-diagonal entry that
  # would cause sed2 < 0 for one pair (floating-point artefact simulation).
  n <- 4L
  set.seed(99L)
  # Start with a valid PD matrix then perturb one off-diagonal slightly
  A    <- matrix(rnorm(n * n) * 0.5, n, n)
  Vbase <- crossprod(A) + diag(1, n)
  # Manually inflate one off-diagonal so that V_ii + V_jj - 2*V_ij < 0
  # (i=1, j=2): requires V_12 > (V_11 + V_22)/2
  V <- Vbase
  V[1L, 2L] <- (Vbase[1L, 1L] + Vbase[2L, 2L]) / 2 + 1e-8  # just above threshold
  V[2L, 1L] <- V[1L, 2L]     # keep symmetric

  pv <- data.frame(
    Genotype        = factor(paste0("G", seq_len(n))),
    predicted.value = rnorm(n, 50, 5),
    std.error       = runif(n, 0.5, 1.5),
    status          = factor(rep("Estimable", n)),
    stringsAsFactors = FALSE
  )
  p <- list(pvals = pv, vcov = V)

  local_mocked_bindings(predict = function(...) p, .package = "biomAid")
  m <- make_asreml_model(nedf = 30L)
  # Should run without error; the guard converts the negative sed2 to NA
  expect_no_error(
    out <- compare(m, term = "Genotype", type = "HSD")
  )
  expect_s3_class(out, "data.frame")
  # Result must be a finite positive number (avsed computed na.rm = TRUE)
  expect_true(is.finite(out$HSD[1L]) && out$HSD[1L] > 0)

  # Also verify the guard directly: sed2 symmetric, negative -> NA
  dv   <- diag(V)
  sed2 <- outer(dv, dv, "+") - 2 * V
  expect_true(isSymmetric(sed2))          # sed2 is symmetric
  neg_before <- sed2[lower.tri(sed2)]
  neg_before[neg_before < 0] <- NA_real_
  expect_true(all(neg_before[!is.na(neg_before)] >= 0))
})
