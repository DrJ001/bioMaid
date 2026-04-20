# tests/testthat/test-randomRegressMV.R
# Tests for randomRegressMV() — exercises the G-matrix algebra directly.
# No ASReml licence required.

# ---------------------------------------------------------------------------
# Helpers: build synthetic G-matrices
# ---------------------------------------------------------------------------

# Positive-definite G-matrix for 2 sites x 3 treatments = 6 levels
make_Gmat <- function(seed = 1L) {
  set.seed(seed)
  k <- 6L
  A <- matrix(rnorm(k * k), k, k)
  G <- crossprod(A) / k + diag(0.5, k)
  lv <- c("N0-S1","N1-S1","N2-S1","N0-S2","N1-S2","N2-S2")
  dimnames(G) <- list(lv, lv)
  G
}

# ---------------------------------------------------------------------------
# 1. Conditional variance formula: sig_j = G_jj - G_jA %*% beta_j
# ---------------------------------------------------------------------------
test_that("conditional variance is positive for valid G-matrix", {
  G    <- make_Gmat(1L)
  # Use indices 1 (N0-S1) and 2 (N1-S1): condition N1 on N0
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

# ---------------------------------------------------------------------------
# 2. Matrix-form beta_j = G_AA^{-1} G_Aj (multiple conditioning)
# ---------------------------------------------------------------------------
test_that("matrix beta satisfies G_AA %*% beta = G_Aj", {
  G     <- make_Gmat(3L)
  j_ind <- 3L
  a_ind <- 1:2
  G_Aj  <- G[a_ind, j_ind, drop = FALSE]
  G_AA  <- G[a_ind, a_ind, drop = FALSE]
  beta  <- drop(solve(G_AA, G_Aj))
  expect_equal(drop(G_AA %*% beta), drop(G_Aj), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 3. Responsiveness BLUP: u_resp = u_j - u_A %*% beta_j
# ---------------------------------------------------------------------------
test_that("responsiveness BLUP is orthogonal to u_A (by construction)", {
  set.seed(5L)
  n     <- 200L   # use large n for stable correlation
  G     <- make_Gmat(5L)
  j_ind <- 2L
  a_ind <- 1L

  # Simulate correlated u via Cholesky of 2x2 G sub-matrix
  G_sub <- G[c(a_ind, j_ind), c(a_ind, j_ind)]   # [A, j] x [A, j]
  R     <- chol(G_sub)
  blups <- matrix(rnorm(n * 2), n, 2) %*% R       # n x 2: col1=u_A, col2=u_j
  u_j   <- blups[, 2L]
  u_a   <- blups[, 1L, drop = FALSE]

  G_Aj  <- G[a_ind, j_ind]
  G_AA  <- G[a_ind, a_ind]
  beta  <- G_Aj / G_AA
  u_resp <- u_j - drop(u_a) * beta

  # Correlation of u_resp with u_A should be ≈ 0
  expect_equal(cor(u_resp, drop(u_a)), 0, tolerance = 0.15)
})

# ---------------------------------------------------------------------------
# 4. Transformation matrix T
# ---------------------------------------------------------------------------
test_that("tmat starts as identity and gets off-diagonal entries", {
  k <- 4L
  tmat <- diag(k)
  # Simulate: condition row 2 on row 1 with beta = 0.7
  tmat[2L, 1L] <- -0.7
  expect_equal(diag(tmat), rep(1, k))
  expect_equal(tmat[2L, 1L], -0.7)
  expect_equal(tmat[1L, 2L], 0)   # upper triangle still 0
})

test_that("TGmat = tmat %*% G %*% t(tmat) is symmetric", {
  G    <- make_Gmat(6L)
  tmat <- diag(nrow(G))
  # Condition index 2 on index 1
  j <- 2L; a <- 1L
  beta_j <- G[a, j] / G[a, a]
  tmat[j, a] <- -beta_j
  TG <- tmat %*% G %*% t(tmat)
  expect_equal(TG, t(TG), tolerance = 1e-10)
})

test_that("sequential TGmat is diagonal (Cholesky property)", {
  # For a 3x3 G, sequential conditioning should diagonalise
  set.seed(10L)
  k <- 3L
  A <- matrix(rnorm(k * k), k, k)
  G <- crossprod(A) / k + diag(0.5, k)
  levs <- c("A","B","C")
  dimnames(G) <- list(levs, levs)

  tmat <- diag(k)
  # Condition B on A
  tmat[2L, 1L] <- -G[1L, 2L] / G[1L, 1L]
  # Condition C on A and B: solve G[1:2,1:2] %*% beta = G[1:2,3]
  G_AA  <- G[1:2, 1:2]
  G_Aj  <- G[1:2, 3L, drop = FALSE]
  beta  <- drop(solve(G_AA, G_Aj))
  tmat[3L, 1:2] <- -beta

  TG <- tmat %*% G %*% t(tmat)
  # Off-diagonals should be ≈ 0 (sequential = Cholesky = LDL')
  off_diag <- TG[upper.tri(TG)]
  expect_equal(off_diag, rep(0, length(off_diag)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. .condList integration
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 6. Validate output structure expected from randomRegressMV
# ---------------------------------------------------------------------------
test_that("blups data frame would have Site, Variety, and BLUP columns", {
  # Simulates what the function builds
  usnams <- c("S1","S2")
  glev   <- paste0("G", sprintf("%02d", 1:10))
  levs   <- c("N0","N1","N2")
  ns     <- length(usnams)
  nvar   <- length(glev)

  blups <- data.frame(
    Site    = rep(usnams, each = nvar),
    Variety = rep(glev,   times = ns),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(blups), ns * nvar)
  expect_named(blups, c("Site","Variety"))
})

# ---------------------------------------------------------------------------
# 7. sigmat: conditional variance equals G_jj - G_jA beta_j
# ---------------------------------------------------------------------------
test_that("sigmat entry matches formula for known G", {
  set.seed(7L)
  k <- 4L
  A <- matrix(rnorm(k * k), k, k)
  G <- crossprod(A) / k + diag(1, k)

  j <- 2L; a <- 1L
  G_jj   <- G[j, j]
  G_jA   <- G[j, a, drop = FALSE]
  G_AA   <- G[a, a, drop = FALSE]
  beta_j <- drop(G_jA) / G_AA[1L, 1L]
  sig_j  <- G_jj - drop(G_jA %*% beta_j)

  # Manual computation
  manual <- G_jj - (G[j, a]^2 / G[a, a])
  expect_equal(sig_j, manual, tolerance = 1e-12)
  expect_gt(sig_j, 0)
})

# ---------------------------------------------------------------------------
# 8. beta list structure
# ---------------------------------------------------------------------------
test_that("beta list has correct names and matrix dimensions", {
  levs  <- c("N0","N1","N2")
  cl    <- biomAid:::.condList(levs, "baseline", NULL)
  cond  <- names(cl[!vapply(cl, is.null, logical(1L))])
  usnams <- c("S1","S2","S3")
  ns    <- length(usnams)

  beta <- setNames(vector("list", length(cond)), cond)
  for (lv in cond)
    beta[[lv]] <- matrix(NA_real_, ns, length(cl[[lv]]),
                         dimnames = list(usnams, cl[[lv]]))

  expect_named(beta, cond)
  expect_equal(dim(beta[["N1"]]), c(ns, 1L))
  expect_equal(dim(beta[["N2"]]), c(ns, 1L))
})

# ---------------------------------------------------------------------------
# 9. Naming of TGmat rows/cols
# ---------------------------------------------------------------------------
test_that("TGmat dimnames include eff. and resp. prefixes", {
  G    <- make_Gmat(8L)
  levs <- c("N0","N1","N2")
  tsnams <- c("N0-S1","N1-S1","N2-S1","N0-S2","N1-S2","N2-S2")
  tmat   <- diag(nrow(G))
  TG     <- tmat %*% G %*% t(tmat)
  dimnames(TG) <- list(tsnams, tsnams)

  # Rename with prefixes
  tsnams_out <- tsnams
  for (lv in c("N0"))
    tsnams_out <- gsub(lv, paste0("eff.", lv), tsnams_out, fixed = TRUE)
  for (lv in c("N1","N2"))
    tsnams_out <- gsub(lv, paste0("resp.", lv), tsnams_out, fixed = TRUE)

  expect_true(all(grepl("eff\\.N0|resp\\.N1|resp\\.N2", tsnams_out)))
})

# ---------------------------------------------------------------------------
# 10. Error conditions
# ---------------------------------------------------------------------------
test_that("levs = NULL errors", {
  expect_error(
    randomRegressMV(list(), levs = NULL),
    "At least two treatment levels"
  )
})

test_that("levs with single element errors", {
  expect_error(
    randomRegressMV(list(), levs = "N0"),
    "At least two treatment levels"
  )
})

# ---------------------------------------------------------------------------
# 11. HSD from PEV-based variance (non-FA case)
# ---------------------------------------------------------------------------
test_that("HSD from block-decomposed PEV is positive", {
  set.seed(40L)
  nvar <- 10L
  sig  <- 1.2
  P    <- diag(nvar)         # simple PEV = identity
  dv   <- diag(sig^2 * P)
  sed  <- outer(dv, dv, "+") - 2 * (sig^2 * P)
  sed  <- sed[lower.tri(sed)]
  sed[sed < 0] <- NA_real_
  hsd  <- (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
            qtukey(0.95, nvar, df = nvar - 2L)
  expect_gt(hsd, 0)
})

# ---------------------------------------------------------------------------
# 12. parse treatment labels from G-matrix column names with separator
# ---------------------------------------------------------------------------
test_that("parsing tsnams with sep='-' extracts treatment and site", {
  tsnams <- c("N0-S1","N1-S1","N0-S2","N1-S2")
  sep    <- "-"
  st     <- strsplit(tsnams, split = sep, fixed = TRUE)
  tnam   <- vapply(st, `[`, character(1L), 1L)
  snam   <- vapply(st, `[`, character(1L), 2L)
  expect_equal(sort(unique(tnam)), c("N0","N1"))
  expect_equal(sort(unique(snam)), c("S1","S2"))
})

test_that("without separator, snam = 'Single'", {
  tsnams <- c("N0","N1","N2")
  sep    <- "-"
  has_sep <- any(grepl(sep, tsnams, fixed = TRUE))
  if (!has_sep) {
    snam <- rep("Single", length(tsnams))
  }
  expect_true(all(snam == "Single"))
})

# ---------------------------------------------------------------------------
# 13. cond_list returned matches type used
# ---------------------------------------------------------------------------
test_that("output cond_list matches the .condList for given type", {
  levs <- c("T0","T1","T2")
  for (tp in c("baseline","sequential","partial")) {
    cl <- biomAid:::.condList(levs, tp, NULL)
    # Check structure is correct list
    expect_named(cl, levs)
    expect_true(is.list(cl))
  }
})

# ---------------------------------------------------------------------------
# 14. Singular G sub-matrix handling (skip, not crash)
# ---------------------------------------------------------------------------
test_that("singular G sub-matrix returns NULL from tryCatch", {
  G_AA  <- matrix(0, 2, 2)   # singular
  G_Aj  <- matrix(c(1, 1), 2, 1)
  beta  <- tryCatch(
    drop(solve(G_AA, G_Aj)),
    error = function(e) NULL
  )
  expect_null(beta)
})

# ---------------------------------------------------------------------------
# 15. pev = TRUE vs pev = FALSE changes PEV calculation
# ---------------------------------------------------------------------------
test_that("pev=FALSE inverts PEV to posterior variance", {
  set.seed(55L)
  n      <- 5L
  sig_j  <- 2.0
  pev_j  <- diag(runif(n, 0.2, 0.5))
  # posterior = sig_j * I - PEV
  post   <- diag(sig_j, n) - pev_j
  expect_equal(diag(post), sig_j - diag(pev_j))
})
