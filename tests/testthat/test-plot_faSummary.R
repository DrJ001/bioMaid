# tests/testthat/test-plot_faSummary.R
# Tests for plot_faSummary() in the biomAid package.
#
# Structure
# =========
#  A. Input validation      – argument checking.
#  B. Rejected / moved types.
#  C. Precondition checks   – k >= 2, VAF data, BLUPs, term selection.
#  D. Return values         – ggplot / list($plot, $data).
#  E. Data structure        – columns and dimensions of $data.
#  F. Internal helpers.
#
# A synthetic faSummary object is used throughout so the file needs no ASReml
# licence.  It is built with the real internal helpers so it stays consistent
# with production output.

# ===========================================================================
# SHARED FIXTURE — 4 environments, 3 varieties, k = 2
# ===========================================================================

.pfs_envs <- c("E1", "E2", "E3", "E4")
.pfs_vars <- c("V1", "V2", "V3")
.pfs_psi  <- c(0.10, 0.20, 0.15, 0.25)
.pfs_lam  <- matrix(c(1.0, 0.8, 1.2, 0.6,
                      0.3, -0.5, 0.4, -0.7), nrow = 4L, ncol = 2L)

.mk_fasummary <- function(lam = .pfs_lam, psi = .pfs_psi,
                          envs = .pfs_envs, vars = .pfs_vars,
                          term = "fa(Site, 2):Variety",
                          with_blups = TRUE) {
  s <- .fa_structure(c(psi, as.vector(lam)), envs)
  v <- .fa_vaf(s$loads, s$spec_var, envs, "Site")
  k <- s$k

  g <- list(
    Gmat = s$Gmat, Cmat = s$Cmat, loads = s$loads, loads_cor = s$loads_cor,
    spec_var = s$spec_var, vaf_env = v$vaf_env, vaf_summary = v$vaf_summary,
    vaf_total = v$vaf_total, k = k, env = envs,
    outer = "Site", inner = "Variety", inner_fun = "id"
  )

  b <- NULL
  if (with_blups) {
    set.seed(11L)
    scores <- data.frame(
      Site    = rep(paste0("Comp", seq_len(k)), each = length(vars)),
      Variety = rep(vars, times = k),
      blup    = stats::rnorm(length(vars) * k),
      stringsAsFactors = FALSE
    )
    scores$blupr <- scores$blup
    score_mat <- matrix(scores$blupr, ncol = k)
    blups <- data.frame(
      Site    = rep(envs, each  = length(vars)),
      Variety = rep(vars, times = length(envs)),
      stringsAsFactors = FALSE
    )
    blups$blup    <- as.vector(score_mat %*% t(s$loads)) + stats::rnorm(nrow(blups), 0, 0.05)
    blups$regblup <- as.vector(score_mat %*% t(s$loads))
    blups$pres    <- 3L
    b <- list(blups = blups, scores = scores)
  }

  structure(
    list(
      gammas = setNames(list(g), term),
      blups  = if (with_blups) setNames(list(b), term) else NULL,
      terms  = term,
      call   = quote(faSummary(model))
    ),
    class = "faSummary"
  )
}

.pfs_res <- .mk_fasummary()

# k = 1 variant
.pfs_res_k1 <- .mk_fasummary(lam = .pfs_lam[, 1L, drop = FALSE],
                             term = "fa(Site, 1):Variety")

# Two-term object with a combined "total" structure carrying no loadings
.mk_fasummary_total <- function() {
  base <- .mk_fasummary()
  tm   <- base$terms[1L]
  g    <- base$gammas[[tm]]
  base$gammas[["Site:Variety-total"]] <- list(
    Gmat = 2 * g$Gmat, Cmat = g$Cmat, env = g$env,
    outer = "Site", inner = "Variety", inner_fun = "total"
  )
  base$blups[["Site:Variety-total"]] <-
    list(blups = base$blups[[tm]]$blups, scores = NULL)
  base$terms <- c(tm, "Site:Variety-total")
  base
}

.pfs_res_tot <- .mk_fasummary_total()

.all_types <- c("VAF", "heatmap", "loadings", "regress", "added")

# ===========================================================================
# SECTION A – Input validation
# ===========================================================================

test_that("A: non-faSummary input is rejected", {
  expect_error(plot_faSummary(data.frame(x = 1)), "must be an object of class")
  expect_error(plot_faSummary(list(a = 1)), "must be an object of class")
})

test_that("A: unknown type is rejected", {
  expect_error(plot_faSummary(.pfs_res, type = "nope"), "should be one of")
})

test_that("A: unknown order is rejected", {
  expect_error(plot_faSummary(.pfs_res, type = "heatmap", order = "nope"),
               "should be one of")
})

test_that("A: non-theme 'theme' is rejected", {
  expect_error(plot_faSummary(.pfs_res, theme = "bw"),
               "must be a ggplot2 theme object")
})

test_that("A: non-logical 'return_data' is rejected", {
  expect_error(plot_faSummary(.pfs_res, return_data = "yes"),
               "single logical value")
  expect_error(plot_faSummary(.pfs_res, return_data = c(TRUE, FALSE)),
               "single logical value")
})

test_that("A: out-of-range 'tol' is rejected", {
  expect_error(plot_faSummary(.pfs_res, type = "loadings", tol = 1.5),
               "between 0 and 1")
  expect_error(plot_faSummary(.pfs_res, type = "loadings", tol = -0.1),
               "between 0 and 1")
  expect_error(plot_faSummary(.pfs_res, type = "loadings", tol = NA),
               "between 0 and 1")
})

test_that("A: unknown term is rejected with available terms listed", {
  expect_error(plot_faSummary(.pfs_res, term = "fa(Nope, 2):Variety"),
               "not found in 'res'")
})

test_that("A: unknown variety is rejected", {
  expect_error(plot_faSummary(.pfs_res, type = "regress",
                              varieties = c("V1", "ZZZ")),
               "variety not found")
})

# ===========================================================================
# SECTION B – Types that live elsewhere
# ===========================================================================

test_that("B: fastIC-only types are rejected", {
  for (typ in c("fast", "biplot", "CVE", "iclass", "OP.pairs", "OP.variety")) {
    expect_error(plot_faSummary(.pfs_res, type = typ), "should be one of",
                 label = typ)
  }
})

# ===========================================================================
# SECTION C – Precondition checks
# ===========================================================================

test_that("C: 'loadings' requires k >= 2", {
  expect_error(plot_faSummary(.pfs_res_k1, type = "loadings"),
               "requires k >= 2", fixed = TRUE)
})

test_that("C: 'VAF' errors when vaf_env is absent", {
  res <- .pfs_res
  res$gammas[[1L]]$vaf_env <- NULL
  expect_error(plot_faSummary(res, type = "VAF"), "VAF data not found")
})

test_that("C: regression plots error when BLUPs are absent", {
  res <- .mk_fasummary(with_blups = FALSE)
  expect_error(plot_faSummary(res, type = "regress"), "not available")
  expect_error(plot_faSummary(res, type = "added"), "not available")
})

test_that("C: the combined total term is rejected for loading-based types", {
  for (typ in c("VAF", "loadings", "regress", "added")) {
    expect_error(
      plot_faSummary(.pfs_res_tot, type = typ, term = "Site:Variety-total"),
      "has no loadings", label = typ
    )
  }
})

test_that("C: the combined total term is allowed for 'heatmap'", {
  p <- plot_faSummary(.pfs_res_tot, type = "heatmap",
                      term = "Site:Variety-total")
  expect_s3_class(p, "ggplot")
})

test_that("C: default term skips terms without loadings", {
  sel <- .pfs_term(.pfs_res_tot, NULL, need_loads = TRUE)
  expect_equal(sel$term, "fa(Site, 2):Variety")
})

# ===========================================================================
# SECTION D – Return values
# ===========================================================================

test_that("D: every type returns a ggplot", {
  for (typ in .all_types) {
    p <- plot_faSummary(.pfs_res, type = typ)
    expect_s3_class(p, "ggplot")
  }
})

test_that("D: return_data = TRUE returns a list with $plot and $data", {
  for (typ in .all_types) {
    out <- plot_faSummary(.pfs_res, type = typ, return_data = TRUE)
    expect_type(out, "list")
    expect_true("plot" %in% names(out), label = paste(typ, "has $plot"))
    expect_true("data" %in% names(out), label = paste(typ, "has $data"))
    expect_s3_class(out$plot, "ggplot")
  }
})

test_that("D: all heatmap orderings return a ggplot", {
  for (o in c("asis", "loading", "cluster")) {
    p <- plot_faSummary(.pfs_res, type = "heatmap", order = o)
    expect_s3_class(p, "ggplot")
  }
})

test_that("D: k = 1 works for VAF, heatmap and the regression plots", {
  for (typ in c("VAF", "heatmap", "regress", "added")) {
    p <- plot_faSummary(.pfs_res_k1, type = typ)
    expect_s3_class(p, "ggplot")
  }
})

test_that("D: explicit varieties are honoured", {
  p <- plot_faSummary(.pfs_res, type = "added", varieties = c("V1", "V3"))
  expect_s3_class(p, "ggplot")
})

# ===========================================================================
# SECTION E – Data structure
# ===========================================================================

test_that("E: 'VAF' $data is a list with $env and $summary", {
  out <- plot_faSummary(.pfs_res, type = "VAF", return_data = TRUE)
  expect_type(out$data, "list")
  expect_named(out$data, c("env", "summary"))
  expect_equal(nrow(out$data$env), length(.pfs_envs))
  expect_equal(nrow(out$data$summary), 3L)   # k = 2 factors + Specific
})

test_that("E: 'VAF' proportions sum to one per environment", {
  out <- plot_faSummary(.pfs_res, type = "VAF", return_data = TRUE)
  env <- out$data$env
  expect_equal(rowSums(env[, c("Factor1", "Factor2", "Specific")]),
               rep(1, length(.pfs_envs)), tolerance = 1e-10)
})

test_that("E: 'heatmap' $data has t^2 rows and a corr column", {
  out <- plot_faSummary(.pfs_res, type = "heatmap", return_data = TRUE)
  expect_equal(nrow(out$data), length(.pfs_envs)^2)
  expect_true(all(c("env_row", "env_col", "corr") %in% names(out$data)))
  expect_true(all(abs(out$data$corr) <= 1 + 1e-12))
})

test_that("E: 'loadings' $data has one row per environment per factor pair", {
  out <- plot_faSummary(.pfs_res, type = "loadings", return_data = TRUE)
  n_pairs <- choose(2L, 2L)
  expect_equal(nrow(out$data), length(.pfs_envs) * n_pairs)
  expect_true(all(c("pair", "lab", "xend", "yend", "radius",
                    "within_tol") %in% names(out$data)))
})

test_that("E: 'loadings' radii never exceed one", {
  out <- plot_faSummary(.pfs_res, type = "loadings", return_data = TRUE)
  expect_true(all(out$data$radius <= 1 + 1e-10))
})

test_that("E: 'tol' controls the within_tol flag", {
  lo <- plot_faSummary(.pfs_res, type = "loadings", tol = 0,
                       return_data = TRUE)$data
  hi <- plot_faSummary(.pfs_res, type = "loadings", tol = 1,
                       return_data = TRUE)$data
  expect_true(all(lo$within_tol))
  expect_true(!any(hi$within_tol))
})

test_that("E: 'regress' $data has variety x factor x environment rows", {
  out <- plot_faSummary(.pfs_res, type = "regress", return_data = TRUE)
  nv  <- length(unique(out$data$variety))
  expect_equal(nrow(out$data), nv * 2L * length(.pfs_envs))
  expect_true(all(c("variety", "facet_lab", "env", "loading",
                    "blup_adj", "slope") %in% names(out$data)))
})

test_that("E: 'added' adjusts the BLUPs but 'regress' does not", {
  reg <- plot_faSummary(.pfs_res, type = "regress", return_data = TRUE)$data
  add <- plot_faSummary(.pfs_res, type = "added",   return_data = TRUE)$data
  expect_equal(nrow(reg), nrow(add))
  expect_false(isTRUE(all.equal(reg$blup_adj, add$blup_adj)))
})

test_that("E: 'added' with k = 1 leaves the BLUPs unchanged", {
  reg <- plot_faSummary(.pfs_res_k1, type = "regress", return_data = TRUE)$data
  add <- plot_faSummary(.pfs_res_k1, type = "added",   return_data = TRUE)$data
  expect_equal(reg$blup_adj, add$blup_adj)
})

# ===========================================================================
# SECTION F – Internal helpers
# ===========================================================================

test_that("F: .pfs_env_order 'asis' preserves the matrix order", {
  g <- .pfs_res$gammas[[1L]]
  expect_equal(.pfs_env_order(g, "asis"), .pfs_envs)
})

test_that("F: .pfs_env_order 'loading' sorts by first-factor loading", {
  g   <- .pfs_res$gammas[[1L]]
  ord <- .pfs_env_order(g, "loading")
  expect_equal(ord, rownames(g$loads)[order(-g$loads[, 1L])])
})

test_that("F: .pfs_env_order 'cluster' returns a permutation", {
  g   <- .pfs_res$gammas[[1L]]
  ord <- .pfs_env_order(g, "cluster")
  expect_setequal(ord, .pfs_envs)
  expect_length(ord, length(.pfs_envs))
})

test_that("F: .pfs_score_matrix has varieties x k dimensions", {
  g <- .pfs_res$gammas[[1L]]
  b <- .pfs_res$blups[[1L]]
  sm <- .pfs_score_matrix(g, b)
  expect_equal(dim(sm), c(length(.pfs_vars), 2L))
  expect_equal(rownames(sm), .pfs_vars)
})

test_that("F: .pfs_default_varieties picks factor extremes and honours n", {
  g  <- .pfs_res$gammas[[1L]]
  sm <- .pfs_score_matrix(g, .pfs_res$blups[[1L]])
  picked <- .pfs_default_varieties(sm, 6L)
  expect_true(all(picked %in% .pfs_vars))
  expect_equal(anyDuplicated(picked), 0L)
  expect_length(.pfs_default_varieties(sm, 2L), 2L)
})
