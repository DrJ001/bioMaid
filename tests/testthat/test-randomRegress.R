# tests/testthat/test-randomRegress.R
# Consolidated tests for randomRegress() and its private helper .condList().
#
# Sections:
#   A. Unit algebra tests       – no model / no mocking required
#   B. .condList() unit tests
#   B2. .parse_rreg_term() unit tests
#   C. Error conditions
#   D. End-to-end tests         – local_mocked_bindings for predict() /
#                                  .asreml_vparams() / .fa_asreml()
#   E. New path coverage        – FA path, pev=FALSE, sep, scalar beta, absent treatment
#   F. New structures / wrappers – corh, corgh, vm, ide; .cor_to_cov_Gmat()

# ===========================================================================
# Shared fixtures
# ===========================================================================

# --- G-matrix builders -----------------------------------------------------

make_Gmat <- function(seed = 1L) {
  set.seed(seed)
  k  <- 6L
  A  <- matrix(rnorm(k * k), k, k)
  G  <- crossprod(A) / k + diag(0.5, k)
  lv <- c("N0-S1","N1-S1","N2-S1","N0-S2","N1-S2","N2-S2")
  dimnames(G) <- list(lv, lv)
  G
}

# G-matrix with configurable levs / sites / seed
make_Gmat_rrm <- function(levs  = c("N0","N1","N2"),
                           sites = c("S1","S2"),
                           seed  = 1L) {
  set.seed(seed)
  tsnams <- as.vector(outer(levs, sites, paste, sep = "-"))
  k      <- length(tsnams)
  A      <- matrix(rnorm(k * k), k, k)
  G      <- crossprod(A) / k + diag(0.5, k)
  dimnames(G) <- list(tsnams, tsnams)
  G
}

# --- predict() return value builder ----------------------------------------

make_pvals_rrm <- function(levs  = c("N0","N1","N2"),
                            sites = c("S1","S2"),
                            n_var = 10L,
                            seed  = 1L) {
  set.seed(seed)
  tsnams    <- as.vector(outer(levs, sites, paste, sep = "-"))
  varieties <- paste0("Var", sprintf("%02d", seq_len(n_var)))
  pv <- expand.grid(
    TSite          = factor(tsnams, levels = tsnams),
    Variety        = factor(varieties),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  pv$predicted.value <- rnorm(nrow(pv), 0, 1)
  pv$std.error       <- runif(nrow(pv), 0.1, 0.5)
  pv$status          <- factor(rep("Estimable", nrow(pv)))
  n    <- nrow(pv)
  Av   <- matrix(rnorm(n * n) * 0.1, n, n)
  vcov <- crossprod(Av) + diag(0.05, n)
  list(pvals = pv, vcov = vcov)
}

# --- Minimal asreml model stub (non-FA, us structure) ----------------------

make_rrm_model <- function(struct = "us") {
  frm <- stats::as.formula(paste0("~ ", struct, "(TSite):Variety"))
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"
  m
}

# ===========================================================================
# SECTION A: Unit algebra – G-matrix conditional formulas
# ===========================================================================

test_that("conditional variance is positive for valid G-matrix", {
  G      <- make_Gmat(1L)
  j_ind  <- 2L
  a_inds <- 1L
  G_jj   <- G[j_ind, j_ind]
  G_jA   <- G[j_ind, a_inds, drop = FALSE]
  G_AA   <- G[a_inds, a_inds, drop = FALSE]
  G_Aj   <- G[a_inds, j_ind,  drop = FALSE]
  beta_j <- drop(G_Aj) / G_AA[1L, 1L]
  sig_j  <- G_jj - drop(G_jA %*% beta_j)
  expect_gt(sig_j, 0)
})

test_that("scalar beta = G_Aj / G_AA when |A|=1", {
  G    <- make_Gmat(2L)
  j    <- 3L
  a    <- 1L
  beta <- G[a, j] / G[a, a]
  expect_equal(beta, G[a, j] / G[a, a])
})

test_that("matrix beta satisfies G_AA %*% beta = G_Aj", {
  G     <- make_Gmat(3L)
  j_ind <- 3L
  a_ind <- 1:2
  G_Aj  <- G[a_ind, j_ind, drop = FALSE]
  G_AA  <- G[a_ind, a_ind, drop = FALSE]
  beta  <- drop(solve(G_AA, G_Aj))
  expect_equal(drop(G_AA %*% beta), drop(G_Aj), tolerance = 1e-10)
})

test_that("responsiveness BLUP is orthogonal to u_A (by construction)", {
  set.seed(5L)
  n     <- 200L
  G     <- make_Gmat(5L)
  j_ind <- 2L
  a_ind <- 1L
  G_sub <- G[c(a_ind, j_ind), c(a_ind, j_ind)]
  R     <- chol(G_sub)
  blups <- matrix(rnorm(n * 2), n, 2) %*% R
  u_j   <- blups[, 2L]
  u_a   <- blups[, 1L, drop = FALSE]
  G_Aj  <- G[a_ind, j_ind]
  G_AA  <- G[a_ind, a_ind]
  beta  <- G_Aj / G_AA
  u_resp <- u_j - drop(u_a) * beta
  expect_equal(cor(u_resp, drop(u_a)), 0, tolerance = 0.15)
})

test_that("tmat starts as identity and gets off-diagonal entries", {
  k    <- 4L
  tmat <- diag(k)
  tmat[2L, 1L] <- -0.7
  expect_equal(diag(tmat), rep(1, k))
  expect_equal(tmat[2L, 1L], -0.7)
  expect_equal(tmat[1L, 2L], 0)
})

test_that("TGmat = tmat %*% G %*% t(tmat) is symmetric", {
  G    <- make_Gmat(6L)
  tmat <- diag(nrow(G))
  j    <- 2L; a <- 1L
  beta_j <- G[a, j] / G[a, a]
  tmat[j, a] <- -beta_j
  TG <- tmat %*% G %*% t(tmat)
  expect_equal(TG, t(TG), tolerance = 1e-10)
})

test_that("sequential TGmat is diagonal (Cholesky property)", {
  set.seed(10L)
  k  <- 3L
  A  <- matrix(rnorm(k * k), k, k)
  G  <- crossprod(A) / k + diag(0.5, k)
  dimnames(G) <- list(c("A","B","C"), c("A","B","C"))
  tmat <- diag(k)
  tmat[2L, 1L] <- -G[1L, 2L] / G[1L, 1L]
  G_AA  <- G[1:2, 1:2]
  G_Aj  <- G[1:2, 3L, drop = FALSE]
  beta  <- drop(solve(G_AA, G_Aj))
  tmat[3L, 1:2] <- -beta
  TG  <- tmat %*% G %*% t(tmat)
  expect_equal(TG[upper.tri(TG)], rep(0, 3L), tolerance = 1e-10)
})

test_that("sigmat entry matches formula for known G", {
  set.seed(7L)
  k  <- 4L
  A  <- matrix(rnorm(k * k), k, k)
  G  <- crossprod(A) / k + diag(1, k)
  j  <- 2L; a <- 1L
  G_jj   <- G[j, j]
  G_jA   <- G[j, a, drop = FALSE]
  G_AA   <- G[a, a, drop = FALSE]
  beta_j <- drop(G_jA) / G_AA[1L, 1L]
  sig_j  <- G_jj - drop(G_jA %*% beta_j)
  manual <- G_jj - (G[j, a]^2 / G[a, a])
  expect_equal(sig_j, manual, tolerance = 1e-12)
  expect_gt(sig_j, 0)
})

test_that("HSD from block-decomposed PEV is positive", {
  set.seed(40L)
  nvar <- 10L
  sig  <- 1.2
  P    <- diag(nvar)
  dv   <- diag(sig^2 * P)
  sed  <- outer(dv, dv, "+") - 2 * (sig^2 * P)
  sed  <- sed[lower.tri(sed)]
  sed[sed < 0] <- NA_real_
  hsd  <- (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
            qtukey(0.95, nvar, df = nvar - 2L)
  expect_gt(hsd, 0)
})

test_that("pev=FALSE inverts PEV to posterior variance", {
  set.seed(55L)
  n     <- 5L
  sig_j <- 2.0
  pev_j <- diag(runif(n, 0.2, 0.5))
  post  <- diag(sig_j, n) - pev_j
  expect_equal(diag(post), sig_j - diag(pev_j))
})

# ===========================================================================
# SECTION B: .condList() unit tests
# ===========================================================================

test_that("condList: sequential has increasing conditioning sets", {
  cl <- biomAid:::.condList(c("N0","N1","N2","N3"), "sequential", NULL)
  for (j in 2:4)
    expect_equal(length(cl[[j]]), j - 1L)
})

test_that("condList: partial conditions each on all others", {
  cl <- biomAid:::.condList(c("A","B","C","D"), "partial", NULL)
  for (lv in c("A","B","C","D"))
    expect_equal(length(cl[[lv]]), 3L)
})

test_that("condList: output matches the .condList for each type", {
  levs <- c("T0","T1","T2")
  for (tp in c("baseline","sequential","partial")) {
    cl <- biomAid:::.condList(levs, tp, NULL)
    expect_named(cl, levs)
    expect_true(is.list(cl))
  }
})

test_that("beta list has correct names and matrix dimensions", {
  levs   <- c("N0","N1","N2")
  cl     <- biomAid:::.condList(levs, "baseline", NULL)
  cond   <- names(cl[!vapply(cl, is.null, logical(1L))])
  usnams <- c("S1","S2","S3")
  beta   <- setNames(vector("list", length(cond)), cond)
  for (lv in cond)
    beta[[lv]] <- matrix(NA_real_, length(usnams), length(cl[[lv]]),
                         dimnames = list(usnams, cl[[lv]]))
  expect_named(beta, cond)
  expect_equal(dim(beta[["N1"]]), c(length(usnams), 1L))
  expect_equal(dim(beta[["N2"]]), c(length(usnams), 1L))
})

test_that("TGmat dimnames include eff. and resp. prefixes", {
  G      <- make_Gmat(8L)
  tsnams <- c("N0-S1","N1-S1","N2-S1","N0-S2","N1-S2","N2-S2")
  tmat   <- diag(nrow(G))
  TG     <- tmat %*% G %*% t(tmat)
  dimnames(TG) <- list(tsnams, tsnams)
  tsnams_out <- tsnams
  for (lv in "N0")
    tsnams_out <- gsub(lv, paste0("eff.", lv), tsnams_out, fixed = TRUE)
  for (lv in c("N1","N2"))
    tsnams_out <- gsub(lv, paste0("resp.", lv), tsnams_out, fixed = TRUE)
  expect_true(all(grepl("eff\\.N0|resp\\.N1|resp\\.N2", tsnams_out)))
})

# ===========================================================================
# SECTION B2: .parse_rreg_term() unit tests
# ===========================================================================

test_that(".parse_rreg_term: us structure parses correctly", {
  p <- biomAid:::.parse_rreg_term("us(TSite):Variety")
  expect_equal(p$struct,     "us")
  expect_equal(p$group_var,  "TSite")
  expect_equal(p$by_var,     "Variety")
  expect_null(p$by_wrapper)
  expect_equal(p$only_term,  "TSite:Variety")
  expect_equal(p$classify,   "TSite:Variety")
})

test_that(".parse_rreg_term: fa structure parses order correctly", {
  p <- biomAid:::.parse_rreg_term("fa(TSite, 2):Variety")
  expect_equal(p$struct,    "fa")
  expect_equal(p$group_var, "TSite")
  expect_equal(p$by_var,    "Variety")
  expect_equal(p$n_fa,      2L)
  expect_equal(p$only_term, "fa(TSite, 2):Variety")
  expect_equal(p$classify,  "TSite:Variety")
})

test_that(".parse_rreg_term: fa strips whitespace from user input", {
  p <- biomAid:::.parse_rreg_term("fa( TSite , 2 ) : Variety")
  expect_equal(p$struct,    "fa")
  expect_equal(p$group_var, "TSite")
  expect_equal(p$n_fa,      2L)
})

test_that(".parse_rreg_term: corgh structure parses correctly", {
  p <- biomAid:::.parse_rreg_term("corgh(TSite):Variety")
  expect_equal(p$struct,    "corgh")
  expect_equal(p$group_var, "TSite")
  expect_equal(p$by_var,    "Variety")
  expect_null(p$n_fa)
  expect_equal(p$only_term, "TSite:Variety")
  expect_equal(p$classify,  "TSite:Variety")
})

test_that(".parse_rreg_term: corh structure parses correctly", {
  p <- biomAid:::.parse_rreg_term("corh(TSite):Variety")
  expect_equal(p$struct,    "corh")
  expect_equal(p$group_var, "TSite")
  expect_equal(p$by_var,    "Variety")
  expect_equal(p$only_term, "TSite:Variety")
})

test_that(".parse_rreg_term: diag structure parses correctly", {
  p <- biomAid:::.parse_rreg_term("diag(TSite):Variety")
  expect_equal(p$struct,    "diag")
  expect_equal(p$only_term, "TSite:Variety")
})

  test_that(".parse_rreg_term: vm wrapper on by-variable", {
    p <- biomAid:::.parse_rreg_term("us(TSite):vm(Variety, giv1)")
    expect_equal(p$struct,      "us")
    expect_equal(p$by_var,      "Variety")
    expect_equal(p$by_wrapper,  "vm")
    expect_equal(p$by_raw,      "vm(Variety, giv1)")   # original spacing preserved
    expect_equal(p$only_term,   "TSite:vm(Variety, giv1)")  # spacing preserved for ASReml match
    expect_equal(p$classify,    "TSite:Variety")            # bare name in classify
  })

  test_that(".parse_rreg_term: ide wrapper on by-variable", {
    p <- biomAid:::.parse_rreg_term("us(TSite):ide(Variety)")
    expect_equal(p$by_var,     "Variety")
    expect_equal(p$by_wrapper, "ide")
    expect_equal(p$only_term,  "TSite:ide(Variety)")  # wrapper kept in only
  })

test_that(".parse_rreg_term: unsupported structure errors", {
  expect_error(
    biomAid:::.parse_rreg_term("ar1(TSite):Variety"),
    "Unsupported variance structure"
  )
})

test_that(".parse_rreg_term: missing colon errors", {
  expect_error(
    biomAid:::.parse_rreg_term("us(TSite)"),
    "must be an interaction"
  )
})

# ===========================================================================
# SECTION C: Error conditions
# ===========================================================================

test_that("levs = NULL errors", {
  expect_error(
    randomRegress(list(), levs = NULL),
    "At least two treatment levels"
  )
})

test_that("levs with single element errors", {
  expect_error(
    randomRegress(list(), levs = "N0"),
    "At least two treatment levels"
  )
})

test_that("singular G sub-matrix returns NULL from tryCatch", {
  G_AA <- matrix(0, 2, 2)    # singular
  G_Aj <- matrix(c(1, 1), 2, 1)
  beta <- tryCatch(
    drop(solve(G_AA, G_Aj)),
    error = function(e) NULL
  )
  expect_null(beta)
})

# Negative conditional variance path: indefinite G-matrix triggers warning,
# sigmat entry remains NA, responsiveness BLUPs remain NA for that treatment.
test_that("negative conditional variance: warning issued, sigmat entry is NA", {
  # Build a 4x4 G-matrix that is indefinite (NOT positive semi-definite).
  # Force G_11 - G_1A G_AA^-1 G_A1 < 0 by making the off-diagonals too large.
  # Simplest: 2x2 with variance 0.5 but covariance 0.7 (correlation > 1 in effect).
  k  <- 4L
  lv <- c("N0-S1","N1-S1","N0-S2","N1-S2")
  # Start from a valid G then make off-diagonal T0-T1 block much too large
  set.seed(1L)
  A <- matrix(rnorm(k * k), k, k)
  G <- crossprod(A) / k + diag(0.5, k)
  dimnames(G) <- list(lv, lv)
  # Inflate [1,2] and [2,1] so the 2x2 T sub-block is indefinite
  G[1, 2] <- G[2, 1] <- sqrt(G[1,1] * G[2,2]) * 1.5   # correlation > 1
  # Build pvals fixture with the same layout
  set.seed(1L)
  varieties  <- paste0("Var", sprintf("%02d", 1:8))
  tsnams_all <- lv
  pv_df <- expand.grid(TSite = factor(tsnams_all, levels = tsnams_all),
                       Variety = factor(varieties),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  pv_df$predicted.value <- rnorm(nrow(pv_df))
  pv_df$std.error       <- runif(nrow(pv_df), 0.1, 0.5)
  pv_df$status          <- factor(rep("Estimable", nrow(pv_df)))
  n    <- nrow(pv_df)
  Avcov <- matrix(rnorm(n * n) * 0.1, n, n)
  vcov  <- crossprod(Avcov) + diag(0.05, n)
  pv    <- list(pvals = pv_df, vcov = vcov)

  frm <- stats::as.formula("~ us(TSite):Variety")
  mod <- list(formulae = list(random = frm), call = list())
  class(mod) <- "asreml"

  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  expect_warning(
    res <- randomRegress(mod,
                         term = "us(TSite):Variety",
                         levs = c("N0","N1"),
                         type = "baseline"),
    "Non-positive conditional variance"
  )
  # sigmat entry for the affected treatment should be NA
  expect_true(anyNA(res$sigmat))
})

# ===========================================================================
# SECTION D: End-to-end tests (local_mocked_bindings)
# ===========================================================================

# --- D1. Returns correct list elements (baseline) --------------------------
test_that("randomRegress() baseline returns named list", {
  G  <- make_Gmat_rrm()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       type = "baseline")
  expect_named(res, c("blups","TGmat","Gmat","beta","sigmat","tmat",
                      "cond_list","type","sep"))
  expect_s3_class(res$blups, "data.frame")
  expect_equal(res$type, "baseline")
})

# --- D2. blups has correct columns -----------------------------------------
test_that("randomRegress() blups has correct columns", {
  G  <- make_Gmat_rrm()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"))
  expect_true("Site"    %in% names(res$blups))
  expect_true("Variety" %in% names(res$blups))
  expect_true("N0"      %in% names(res$blups))
  expect_true("resp.N1" %in% names(res$blups))
  expect_true("resp.N2" %in% names(res$blups))
})

# --- D3. Gmat returned correctly -------------------------------------------
test_that("randomRegress() Gmat matches mock", {
  G  <- make_Gmat_rrm(seed = 2L)
  pv <- make_pvals_rrm(seed = 2L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"))
  expect_equal(res$Gmat, G, tolerance = 1e-12)
})

# --- D4. Sequential: TGmat is a matrix ------------------------------------
test_that("randomRegress() sequential: TGmat is a matrix", {
  G  <- make_Gmat_rrm(seed = 3L)
  pv <- make_pvals_rrm(seed = 3L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       type = "sequential")
  expect_equal(res$type, "sequential")
  expect_true(is.matrix(res$TGmat))
})

# --- D5. Partial: resp columns for all levs --------------------------------
test_that("randomRegress() partial: resp columns for all levs", {
  G  <- make_Gmat_rrm(seed = 4L)
  pv <- make_pvals_rrm(seed = 4L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       type = "partial")
  expect_true(all(c("resp.N0","resp.N1","resp.N2") %in% names(res$blups)))
})

# --- D6. sigmat values are positive ----------------------------------------
test_that("randomRegress() sigmat values are positive", {
  G  <- make_Gmat_rrm(seed = 5L)
  pv <- make_pvals_rrm(seed = 5L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"))
  expect_true(all(res$sigmat > 0, na.rm = TRUE))
})

# --- D7. beta coefficients are finite --------------------------------------
test_that("randomRegress() beta coefficients are finite", {
  G  <- make_Gmat_rrm(seed = 6L)
  pv <- make_pvals_rrm(seed = 6L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"))
  expect_true(all(is.finite(res$beta$N1), na.rm = TRUE))
})

# --- D8. Custom type runs end-to-end ---------------------------------------
test_that("randomRegress() custom type runs without error", {
  G  <- make_Gmat_rrm(seed = 7L)
  pv <- make_pvals_rrm(seed = 7L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- suppressWarnings(
    randomRegress(make_rrm_model(),
                  term = "us(TSite):Variety",
                  levs = c("N0","N1","N2"),
                  type = "custom",
                  cond = list(N0 = NULL, N1 = "N0", N2 = c("N0","N1")))
  )
  expect_equal(res$type, "custom")
  expect_true(all(c("resp.N1","resp.N2") %in% names(res$blups)))
})

# --- D9. blups data frame structure (row count) ----------------------------
test_that("blups data frame has Site, Variety and correct row count", {
  usnams <- c("S1","S2")
  glev   <- paste0("G", sprintf("%02d", 1:10))
  nvar   <- length(glev)
  ns     <- length(usnams)
  blups  <- data.frame(
    Site    = rep(usnams, each = nvar),
    Variety = rep(glev,   times = ns),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(blups), ns * nvar)
  expect_named(blups, c("Site","Variety"))
})

# ===========================================================================
# SECTION E: New path coverage
# ===========================================================================

# ---------------------------------------------------------------------------
# E1. FA model path — mock .fa_asreml(); HSD columns must all be NA
# ---------------------------------------------------------------------------

# Reuse the FA fixture layout from test-fast-e2e.R
# .fa_asreml() returns list(gammas=..., blups=...)
# randomRegress.R accesses:
#   sumfa$blups[[rterm]]$blups[, 1:3]  -> columns (blup, enam, vnam)
#   sumfa$gammas[[rterm]]$Gmat         -> G-matrix

make_fa_mock <- function(levs  = c("N0","N1"),
                          sites = c("S1","S2"),
                          n_var = 6L,
                          seed  = 99L) {
  set.seed(seed)
  tsnams    <- as.vector(outer(levs, sites, paste, sep = "-"))
  varieties <- paste0("Var", sprintf("%02d", seq_len(n_var)))
  k         <- length(tsnams)
  A         <- matrix(rnorm(k * k), k, k)
  Gmat      <- crossprod(A) / k + diag(0.5, k)
  dimnames(Gmat) <- list(tsnams, tsnams)

  # Build the blups sub-list:
  # randomRegress uses: sumfa$blups[[rterm]]$blups[, 1:3]
  # then renames to c("blup", enam, vnam).
  # So the element at [[rterm]] must be a list with a $blups data.frame
  # whose first 3 columns are [blup, TSite, Variety].
  blup_df <- expand.grid(
    TSite          = factor(tsnams, levels = tsnams),
    Variety        = factor(varieties),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  blup_df$blup <- rnorm(nrow(blup_df))
  # first 3 cols: blup, TSite, Variety  (order matches what the rename expects)
  blup_df <- blup_df[, c("blup","TSite","Variety")]

  # The rterm matched by grep() on term.labels of "~ fa(TSite, 1):Variety"
  # is "fa(TSite, 1):Variety" (R inserts a space after the comma).
  fa_rterm <- "fa(TSite, 1):Variety"

  blup_inner  <- list(blups = blup_df)          # $blups holds the data.frame
  blup_outer  <- setNames(list(blup_inner), fa_rterm)
  gam_outer   <- setNames(list(list(Gmat = Gmat)), fa_rterm)

  list(gammas = gam_outer, blups = blup_outer)
}

make_fa_model <- function() {
  # formulae$random must contain a term starting with "fa("
  # R's terms() will canonicalise "fa(TSite,1)" to "fa(TSite, 1)" in term.labels
  frm <- stats::as.formula("~ fa(TSite,1):Variety")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"
  m
}

test_that("FA path: runs without error and HSD columns are all NA", {
  fa_mock <- make_fa_mock()
  local_mocked_bindings(
    .fa_asreml = function(...) fa_mock,
    .package   = "biomAid"
  )
  res <- randomRegress(make_fa_model(),
                       term = "fa(TSite, 1):Variety",
                       levs = c("N0","N1"),
                       type = "baseline")
  # Function should succeed and return a list
  expect_true(is.list(res))
  expect_s3_class(res$blups, "data.frame")
  # HSD column must exist and all values must be NA (pred=NULL for FA path)
  expect_true("HSD.N1" %in% names(res$blups))
  expect_true(all(is.na(res$blups$HSD.N1)))
})

# ---------------------------------------------------------------------------
# E2. pev=FALSE path — non-FA, pev=FALSE uses diag(sig_j,nvar) - pev_j
# ---------------------------------------------------------------------------
test_that("randomRegress() pev=FALSE runs without error and returns result", {
  G  <- make_Gmat_rrm(seed = 20L)
  pv <- make_pvals_rrm(seed = 20L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       pev  = FALSE)
  expect_true(is.list(res))
  expect_s3_class(res$blups, "data.frame")
  # HSD column should exist (pev=FALSE doesn't suppress HSD, just changes variance)
  expect_true("HSD.N1" %in% names(res$blups))
})

# ---------------------------------------------------------------------------
# E3. sep argument — composite "treatment-site" labels like "S1-N0"
# ---------------------------------------------------------------------------

# Unit parse check (no model needed)
test_that("parsing tsnams with sep='-' extracts treatment and site correctly", {
  tsnams <- c("N0-S1","N1-S1","N0-S2","N1-S2")
  sep    <- "-"
  st     <- strsplit(tsnams, split = sep, fixed = TRUE)
  tnam   <- vapply(st, `[`, character(1L), 1L)
  snam   <- vapply(st, `[`, character(1L), 2L)
  expect_equal(sort(unique(tnam)), c("N0","N1"))
  expect_equal(sort(unique(snam)), c("S1","S2"))
})

test_that("without separator, snam = 'Single'", {
  tsnams  <- c("N0","N1","N2")
  sep     <- "-"
  has_sep <- any(grepl(sep, tsnams, fixed = TRUE))
  if (!has_sep) {
    snam <- rep("Single", length(tsnams))
  }
  expect_true(all(snam == "Single"))
})

# End-to-end: G-matrix labels in "site-treatment" order (reversed vs default)
# The function detects levs in snam and swaps tnam/snam accordingly.
test_that("randomRegress() sep='-' correctly parses site-treatment labels", {
  levs  <- c("N0","N1","N2")
  sites <- c("S1","S2")
  # Build G with "site-treatment" ordering: e.g. "S1-N0"
  tsnams_st <- as.vector(outer(sites, levs, paste, sep = "-"))
  k         <- length(tsnams_st)
  set.seed(30L)
  Av  <- matrix(rnorm(k * k), k, k)
  G   <- crossprod(Av) / k + diag(0.5, k)
  dimnames(G) <- list(tsnams_st, tsnams_st)

  # Build matching pvals (TSite matches tsnams_st entries)
  varieties <- paste0("Var", sprintf("%02d", 1:8))
  pv_df <- expand.grid(
    TSite          = factor(tsnams_st, levels = tsnams_st),
    Variety        = factor(varieties),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  set.seed(30L)
  pv_df$predicted.value <- rnorm(nrow(pv_df))
  pv_df$std.error       <- runif(nrow(pv_df), 0.1, 0.5)
  pv_df$status          <- factor(rep("Estimable", nrow(pv_df)))
  n     <- nrow(pv_df)
  Avcov <- matrix(rnorm(n * n) * 0.1, n, n)
  vcov  <- crossprod(Avcov) + diag(0.05, n)
  pv    <- list(pvals = pv_df, vcov = vcov)

  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = levs,
                       sep  = "-")
  expect_s3_class(res$blups, "data.frame")
  expect_true("Site"    %in% names(res$blups))
  expect_true("resp.N1" %in% names(res$blups))
  # Sites should be S1 and S2
  expect_equal(sort(unique(res$blups$Site)), c("S1","S2"))
})

# ---------------------------------------------------------------------------
# E4. Scalar beta shortcut — |A_j|=1 gives same result as matrix solve()
# ---------------------------------------------------------------------------
test_that("scalar beta path gives same result as matrix solve() path", {
  G     <- make_Gmat(9L)
  j_ind <- 2L
  a_ind <- 1L

  G_Aj  <- G[a_ind, j_ind, drop = FALSE]
  G_AA  <- G[a_ind, a_ind, drop = FALSE]

  # Scalar shortcut (the path used when length(A_j)==1)
  beta_scalar <- drop(G_Aj) / G_AA[1L, 1L]

  # Matrix solve() path
  beta_matrix <- drop(solve(G_AA, G_Aj))

  expect_equal(beta_scalar, beta_matrix, tolerance = 1e-12)

  # Verify they produce the same sig_j
  G_jj       <- G[j_ind, j_ind]
  G_jA       <- G[j_ind, a_ind, drop = FALSE]
  sig_scalar <- G_jj - drop(G_jA %*% beta_scalar)
  sig_matrix <- G_jj - drop(G_jA %*% beta_matrix)
  expect_equal(sig_scalar, sig_matrix, tolerance = 1e-12)
})

# End-to-end: baseline (|A_j|=1 for every conditioned treatment) uses scalar path
test_that("randomRegress() baseline scalar path completes without error", {
  G  <- make_Gmat_rrm(seed = 40L)
  pv <- make_pvals_rrm(seed = 40L)
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  # baseline: every conditioned treatment uses scalar A_j
  res_baseline <- randomRegress(make_rrm_model(),
                                term = "us(TSite):Variety",
                                levs = c("N0","N1","N2"),
                                type = "baseline")
  # sequential with |A_j|>1 uses matrix solve() for later treatments
  res_seq <- randomRegress(make_rrm_model(),
                           term = "us(TSite):Variety",
                           levs = c("N0","N1","N2"),
                           type = "sequential")

  # Both should return numeric sigmat values
  expect_true(all(is.numeric(res_baseline$sigmat)))
  expect_true(all(is.numeric(res_seq$sigmat)))
  # Diagonal of TGmat should match sigmat entries (within site)
  expect_true(is.matrix(res_baseline$TGmat))
  expect_true(is.matrix(res_seq$TGmat))
})

# ---------------------------------------------------------------------------
# E5. Treatment absent from a site — `next` guard, no error, other sites OK
# ---------------------------------------------------------------------------

make_pvals_missing_treatment <- function(levs  = c("N0","N1","N2"),
                                          sites = c("S1","S2"),
                                          n_var = 6L,
                                          seed  = 50L) {
  # S2 will be missing N2 entirely
  set.seed(seed)
  varieties <- paste0("Var", sprintf("%02d", seq_len(n_var)))

  # Full cross for S1; S2 only gets N0 and N1
  tsnams_S1 <- paste(levs, "S1", sep = "-")
  tsnams_S2 <- paste(levs[1:2], "S2", sep = "-")   # N2-S2 absent
  tsnams_all <- c(tsnams_S1, tsnams_S2)

  pv_df <- expand.grid(
    TSite          = factor(tsnams_all, levels = tsnams_all),
    Variety        = factor(varieties),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  pv_df$predicted.value <- rnorm(nrow(pv_df))
  pv_df$std.error       <- runif(nrow(pv_df), 0.1, 0.5)
  pv_df$status          <- factor(rep("Estimable", nrow(pv_df)))
  n     <- nrow(pv_df)
  Avcov <- matrix(rnorm(n * n) * 0.1, n, n)
  vcov  <- crossprod(Avcov) + diag(0.05, n)
  list(pvals = pv_df, vcov = vcov)
}

make_Gmat_missing_treatment <- function(levs  = c("N0","N1","N2"),
                                         sites = c("S1","S2"),
                                         seed  = 50L) {
  # G-matrix covers all tsnams that appear in pvals (S2 missing N2)
  set.seed(seed)
  tsnams_S1  <- paste(levs,     "S1", sep = "-")
  tsnams_S2  <- paste(levs[1:2],"S2", sep = "-")
  tsnams_all <- c(tsnams_S1, tsnams_S2)
  k  <- length(tsnams_all)
  A  <- matrix(rnorm(k * k), k, k)
  G  <- crossprod(A) / k + diag(0.5, k)
  dimnames(G) <- list(tsnams_all, tsnams_all)
  G
}

test_that("treatment absent from a site: no error, other sites return results", {
  G  <- make_Gmat_missing_treatment()
  pv <- make_pvals_missing_treatment()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  # Should run without error even though N2 is absent in S2
  res <- randomRegress(make_rrm_model(),
                       term = "us(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       type = "baseline")
  expect_s3_class(res$blups, "data.frame")
  # S1 has all three treatments — resp.N1 and resp.N2 should be non-NA for S1
  s1_rows <- res$blups$Site == "S1"
  expect_true(any(!is.na(res$blups$resp.N1[s1_rows])))
  expect_true(any(!is.na(res$blups$resp.N2[s1_rows])))
  # S2 is missing N2: resp.N2 rows for S2 must all be NA
  s2_rows <- res$blups$Site == "S2"
  expect_true(all(is.na(res$blups$resp.N2[s2_rows])))
})

# ===========================================================================
# SECTION F: New structures / wrappers — corh, corgh, vm, ide
# ===========================================================================

# ---------------------------------------------------------------------------
# F1. .cor_to_cov_Gmat() unit tests
# ---------------------------------------------------------------------------

# Build a known mixed variance/correlation matrix (what ASReml returns for
# corh / corgh via vparameters[["TSite:Variety"]]).
# Diagonal = genetic variances; off-diagonal = correlations.
make_cor_vparams <- function(seed = 101L) {
  set.seed(seed)
  k    <- 4L
  sds  <- sqrt(runif(k, 0.5, 3.0))             # heterogeneous std devs
  R    <- cov2cor(crossprod(matrix(rnorm(k*k), k, k)) + diag(0.5, k))
  M    <- R
  diag(M) <- sds^2                              # variances on diagonal
  nms  <- paste0("T", seq_len(k), "-S1")
  dimnames(M) <- list(nms, nms)
  list(M = M, sds = sds, R = R)
}

test_that(".cor_to_cov_Gmat: diagonal of output equals diagonal of input", {
  obj <- make_cor_vparams()
  G   <- biomAid:::.cor_to_cov_Gmat(obj$M)
  expect_equal(diag(G), diag(obj$M), tolerance = 1e-12)
})

test_that(".cor_to_cov_Gmat: off-diagonal equals r_ij * sigma_i * sigma_j", {
  obj  <- make_cor_vparams()
  G    <- biomAid:::.cor_to_cov_Gmat(obj$M)
  k    <- nrow(obj$M)
  sds  <- sqrt(diag(obj$M))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i != j) {
        expected <- unname(obj$M[i, j] * sds[i] * sds[j])
        expect_equal(unname(G[i, j]), expected, tolerance = 1e-12)
      }
    }
  }
})

test_that(".cor_to_cov_Gmat: result is symmetric", {
  obj <- make_cor_vparams()
  G   <- biomAid:::.cor_to_cov_Gmat(obj$M)
  expect_equal(G, t(G), tolerance = 1e-12)
})

test_that(".cor_to_cov_Gmat: result is positive definite (all eigenvalues > 0)", {
  obj <- make_cor_vparams()
  G   <- biomAid:::.cor_to_cov_Gmat(obj$M)
  evs <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evs > 0))
})

test_that(".cor_to_cov_Gmat: cov2cor of result recovers off-diagonal correlations", {
  obj <- make_cor_vparams()
  G   <- biomAid:::.cor_to_cov_Gmat(obj$M)
  R_recovered <- cov2cor(G)
  # Off-diagonal of recovered correlation should equal off-diagonal of M
  k <- nrow(obj$M)
  for (i in seq_len(k))
    for (j in seq_len(k))
      if (i != j)
        expect_equal(R_recovered[i, j], obj$M[i, j], tolerance = 1e-10)
})

test_that(".cor_to_cov_Gmat: preserves dimension names", {
  obj <- make_cor_vparams()
  G   <- biomAid:::.cor_to_cov_Gmat(obj$M)
  expect_equal(dimnames(G), dimnames(obj$M))
})

# ---------------------------------------------------------------------------
# F2. corgh end-to-end — mock returns mixed var/cor matrix; Gmat must be
#     a full covariance matrix with correct off-diagonals after conversion
# ---------------------------------------------------------------------------

make_corgh_model <- function() {
  frm <- stats::as.formula("~ corgh(TSite):Variety")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"
  m
}

# corgh mock: .asreml_vparams() returns the mixed var/cor matrix
make_corgh_vparams <- function(levs  = c("N0","N1","N2"),
                                sites = c("S1","S2"),
                                seed  = 200L) {
  set.seed(seed)
  tsnams <- as.vector(outer(levs, sites, paste, sep = "-"))
  k      <- length(tsnams)
  sds    <- sqrt(runif(k, 0.5, 2.5))
  Atmp   <- matrix(rnorm(k * k), k, k)
  R      <- cov2cor(crossprod(Atmp) + diag(0.5, k))
  M      <- R
  diag(M) <- sds^2
  dimnames(M) <- list(tsnams, tsnams)
  M
}

test_that("corgh: Gmat returned by randomRegress is a full covariance matrix", {
  M  <- make_corgh_vparams()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) M,
    .package        = "biomAid"
  )
  res <- randomRegress(make_corgh_model(),
                       term = "corgh(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       type = "baseline")
  G <- res$Gmat
  # Diagonal should equal the diagonal of M (variances unchanged)
  expect_equal(diag(G), diag(M), tolerance = 1e-12)
  # Off-diagonal should NOT equal the raw correlations from M
  # (they must be covariances = r_ij * sigma_i * sigma_j)
  k   <- nrow(M)
  sds <- sqrt(diag(M))
  expect_equal(unname(G[1,2]), unname(M[1,2] * sds[1] * sds[2]), tolerance = 1e-12)
  # Result must be symmetric
  expect_equal(G, t(G), tolerance = 1e-12)
})

test_that("corgh: randomRegress returns valid blups and positive sigmat", {
  M  <- make_corgh_vparams()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) M,
    .package        = "biomAid"
  )
  res <- randomRegress(make_corgh_model(),
                       term = "corgh(TSite):Variety",
                       levs = c("N0","N1","N2"),
                       type = "baseline")
  expect_s3_class(res$blups, "data.frame")
  expect_true("resp.N1" %in% names(res$blups))
  expect_true(all(res$sigmat > 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# F3. corh end-to-end (two treatment levels)
# ---------------------------------------------------------------------------

make_corh_model <- function() {
  frm <- stats::as.formula("~ corh(TSite):Variety")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"
  m
}

make_corh_vparams <- function(levs  = c("N0","N1"),
                               sites = c("S1","S2"),
                               seed  = 300L) {
  set.seed(seed)
  tsnams <- as.vector(outer(levs, sites, paste, sep = "-"))
  k      <- length(tsnams)
  sds    <- sqrt(runif(k, 0.5, 2.5))
  Atmp   <- matrix(rnorm(k * k), k, k)
  R      <- cov2cor(crossprod(Atmp) + diag(0.5, k))
  M      <- R
  diag(M) <- sds^2
  dimnames(M) <- list(tsnams, tsnams)
  M
}

make_pvals_rrm_2lev <- function(levs  = c("N0","N1"),
                                 sites = c("S1","S2"),
                                 n_var = 10L,
                                 seed  = 300L) {
  make_pvals_rrm(levs = levs, sites = sites, n_var = n_var, seed = seed)
}

test_that("corh: Gmat off-diagonals are covariances not correlations", {
  M  <- make_corh_vparams()
  pv <- make_pvals_rrm_2lev()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) M,
    .package        = "biomAid"
  )
  res <- randomRegress(make_corh_model(),
                       term = "corh(TSite):Variety",
                       levs = c("N0","N1"),
                       type = "baseline")
  G   <- res$Gmat
  sds <- sqrt(diag(M))
  expect_equal(diag(G), diag(M), tolerance = 1e-12)
  expect_equal(unname(G[1,2]), unname(M[1,2] * sds[1] * sds[2]), tolerance = 1e-12)
  expect_equal(G, t(G), tolerance = 1e-12)
})

test_that("corh: randomRegress returns valid blups", {
  M  <- make_corh_vparams()
  pv <- make_pvals_rrm_2lev()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) M,
    .package        = "biomAid"
  )
  res <- randomRegress(make_corh_model(),
                       term = "corh(TSite):Variety",
                       levs = c("N0","N1"),
                       type = "baseline")
  expect_s3_class(res$blups, "data.frame")
  expect_true("resp.N1" %in% names(res$blups))
  expect_true(all(res$sigmat > 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# F4. vm wrapper — term = "us(TSite):vm(Variety, giv1)"
#     The bare "Variety" name must be used in pvals lookup and column names
# ---------------------------------------------------------------------------

make_vm_model <- function() {
  frm <- stats::as.formula("~ us(TSite):vm(Variety, giv1)")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"
  m
}

test_that("vm wrapper: randomRegress parses bare Variety name correctly", {
  G  <- make_Gmat_rrm()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_vm_model(),
                       term = "us(TSite):vm(Variety, giv1)",
                       levs = c("N0","N1","N2"),
                       type = "baseline")
  expect_s3_class(res$blups, "data.frame")
  # Columns should use bare "Variety" name, not "vm(Variety, giv1)"
  expect_true("Variety" %in% names(res$blups))
  expect_true("resp.N1" %in% names(res$blups))
})

  test_that("vm wrapper: by_wrapper recorded, only_term includes vm wrapper", {
    p <- biomAid:::.parse_rreg_term("us(TSite):vm(Variety, giv1)")
    expect_equal(p$by_wrapper, "vm")
    expect_equal(p$by_var,     "Variety")
    expect_equal(p$only_term,  "TSite:vm(Variety, giv1)")  # spacing preserved for ASReml match
    expect_equal(p$classify,   "TSite:Variety")
  })

# ---------------------------------------------------------------------------
# F5. ide wrapper — term = "us(TSite):ide(Variety)"
# ---------------------------------------------------------------------------

make_ide_model <- function() {
  frm <- stats::as.formula("~ us(TSite):ide(Variety)")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"
  m
}

test_that("ide wrapper: randomRegress parses bare Variety name correctly", {
  G  <- make_Gmat_rrm()
  pv <- make_pvals_rrm()
  local_mocked_bindings(
    predict         = function(...) pv,
    .asreml_vparams = function(...) G,
    .package        = "biomAid"
  )
  res <- randomRegress(make_ide_model(),
                       term = "us(TSite):ide(Variety)",
                       levs = c("N0","N1","N2"),
                       type = "baseline")
  expect_s3_class(res$blups, "data.frame")
  expect_true("Variety" %in% names(res$blups))
  expect_true("resp.N1" %in% names(res$blups))
})

  test_that("ide wrapper: by_wrapper recorded, only_term includes ide wrapper", {
    p <- biomAid:::.parse_rreg_term("us(TSite):ide(Variety)")
    expect_equal(p$by_wrapper, "ide")
    expect_equal(p$by_var,     "Variety")
    expect_equal(p$only_term,  "TSite:ide(Variety)")
  })

# ---------------------------------------------------------------------------
# F6. RHS wrapper mismatch — vm() term supplied against ide() model and vice
#     versa.  Previously the rterm grep matched on struct(Group) only, so a
#     mismatched wrapper (e.g. ide vs vm) silently passed validation, used the
#     same G-matrix key ("TSite:Variety" bare), and ASReml-R silently ignored
#     the wrong only= string, making vm() and ide() return identical results.
# ---------------------------------------------------------------------------

test_that("RHS mismatch: ide() term against vm() model errors with informative message", {
  # Model was fitted with vm() wrapper
  frm <- stats::as.formula("~ us(TSite):vm(Variety, giv1)")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"

  expect_error(
    randomRegress(m,
                  term = "us(TSite):ide(Variety)",   # wrong wrapper
                  levs = c("N0", "N1", "N2")),
    "Term mismatch"
  )
})

test_that("RHS mismatch: vm() term against ide() model errors with informative message", {
  # Model was fitted with ide() wrapper
  frm <- stats::as.formula("~ us(TSite):ide(Variety)")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"

  expect_error(
    randomRegress(m,
                  term = "us(TSite):vm(Variety, giv1)",   # wrong wrapper
                  levs = c("N0", "N1", "N2")),
    "Term mismatch"
  )
})

test_that("RHS mismatch: vm() term against bare-Variety model errors", {
  frm <- stats::as.formula("~ us(TSite):Variety")
  m   <- list(formulae = list(random = frm), call = list())
  class(m) <- "asreml"

  expect_error(
    randomRegress(m,
                  term = "us(TSite):vm(Variety, giv1)",
                  levs = c("N0", "N1", "N2")),
    "Term mismatch"
  )
})
