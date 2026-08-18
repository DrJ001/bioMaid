# tests/testthat/test-faSummary.R
# Tests for faSummary() in the biomAid package.
#
# Structure
# =========
#  A. Internal maths        – pure-R unit tests, no ASReml licence needed.
#  B. Golden-value tests    – faSummary() vs the real ASExtras4::fa.asreml().
#  C. API / structure tests – term filtering, blups=FALSE, combine.ide, print.
#  D. Error handling        – validation messages.
#
# Sections B and C require a licensed ASReml-R installation and the fixtures
# built by scratch/make_fa_fixtures.R.  They skip cleanly when unavailable.
# NOTE: reading the fixtures loads the asreml namespace, so the skip guard
# must come before readRDS().

# ---------------------------------------------------------------------------
# Fixture access
# ---------------------------------------------------------------------------
.fix_path <- function(name) test_path("fixtures", paste0(name, ".rds"))

.skip_unless_fixture <- function(name) {
  skip_if_not_installed("asreml")
  skip_if_not(file.exists(.fix_path(name)),
              paste0("fixture ", name, ".rds not built"))
}

.load_fix <- function(name) readRDS(.fix_path(name))

# ===========================================================================
# SECTION A – Internal maths (pure R, no ASReml)
# ===========================================================================

# 3 environments, 2 factors.  matrix(vpar, nrow = 3) fills column-major, so
# the layout is c(psi, loading column 1, loading column 2).
.psi_a  <- c(0.10, 0.20, 0.05)
.lam_a  <- matrix(c(1.0, 0.8, 1.2,
                    0.3, -0.5, 0.4), nrow = 3L, ncol = 2L)
.vpar_a <- c(.psi_a, as.vector(.lam_a))
.envs_a <- c("E1", "E2", "E3")

test_that("A: .fa_structure returns correct k and dimensions", {
  s <- .fa_structure(.vpar_a, .envs_a)
  expect_equal(s$k, 2L)
  expect_equal(dim(s$loads), c(3L, 2L))
  expect_equal(dim(s$Gmat), c(3L, 3L))
  expect_equal(rownames(s$loads), .envs_a)
  expect_equal(colnames(s$loads), c("fac_1", "fac_2"))
})

test_that("A: rotation leaves the covariance structure invariant", {
  s <- .fa_structure(.vpar_a, .envs_a)
  # Lambda %*% t(Lambda) is invariant to any orthogonal rotation
  expect_equal(unname(s$loads %*% t(s$loads)),
               unname(.lam_a %*% t(.lam_a)), tolerance = 1e-12)
})

test_that("A: Gmat = loads %*% t(loads) + diag(spec_var)", {
  s <- .fa_structure(.vpar_a, .envs_a)
  expect_equal(unname(s$Gmat),
               unname(s$loads %*% t(s$loads) + diag(.psi_a)),
               tolerance = 1e-12)
})

test_that("A: Cmat is the correlation matrix of Gmat", {
  s <- .fa_structure(.vpar_a, .envs_a)
  expect_equal(unname(s$Cmat), unname(stats::cov2cor(s$Gmat)),
               tolerance = 1e-12)
  expect_equal(unname(diag(s$Cmat)), rep(1, 3L), tolerance = 1e-12)
})

test_that("A: rotation uses the ASExtras4 convention (negated SVD)", {
  s <- .fa_structure(.vpar_a, .envs_a)
  expect_equal(unname(s$loads), unname(-.lam_a %*% svd(.lam_a)$v),
               tolerance = 1e-12)
})

test_that("A: k = 1 leaves loadings unrotated and rotation NULL", {
  vpar <- c(.psi_a, c(1.0, 0.8, 1.2))
  s <- .fa_structure(vpar, .envs_a)
  expect_equal(s$k, 1L)
  expect_null(s$rotation)
  expect_equal(unname(as.vector(s$loads)), c(1.0, 0.8, 1.2))
})

test_that("A: loads_cor rows satisfy sum of squares = 1 - psi/diag(G)", {
  s <- .fa_structure(.vpar_a, .envs_a)
  expect_equal(unname(rowSums(s$loads_cor^2)),
               unname(1 - .psi_a / diag(s$Gmat)), tolerance = 1e-12)
})

test_that("A: .fa_vaf rows sum to one", {
  s <- .fa_structure(.vpar_a, .envs_a)
  v <- .fa_vaf(s$loads, s$spec_var, .envs_a, "Site")
  props <- v$vaf_env[, c("Factor1", "Factor2", "Specific")]
  expect_equal(unname(rowSums(props)), rep(1, 3L), tolerance = 1e-12)
})

test_that("A: .fa_vaf summary proportions sum to one and cumulate", {
  s <- .fa_structure(.vpar_a, .envs_a)
  v <- .fa_vaf(s$loads, s$spec_var, .envs_a, "Site")
  expect_equal(sum(v$vaf_summary$pct_var), 1, tolerance = 1e-12)
  expect_equal(v$vaf_summary$cum_pct, cumsum(v$vaf_summary$pct_var))
  expect_equal(v$vaf_summary$factor, c("Factor 1", "Factor 2", "Specific"))
})

test_that("A: .fa_vaf column layout and names are stable", {
  s <- .fa_structure(.vpar_a, .envs_a)
  v <- .fa_vaf(s$loads, s$spec_var, .envs_a, "Site")
  expect_equal(names(v$vaf_env),
               c("Site", "Factor1", "Factor2", "Specific", "total_var"))
  expect_equal(v$vaf_env$Site, .envs_a)
})

test_that("A: .fa_vaf clamps negative boundary specific variances", {
  s <- .fa_structure(.vpar_a, .envs_a)
  neg <- c(-1e-9, 0.20, 0.05)
  v <- .fa_vaf(s$loads, neg, .envs_a, "Site")
  expect_true(all(v$vaf_env$Specific >= 0))
  expect_equal(v$vaf_env$Specific[1L], 0, tolerance = 1e-12)
})

test_that("A: vaf_total is the percentage explained by all factors", {
  s <- .fa_structure(.vpar_a, .envs_a)
  v <- .fa_vaf(s$loads, s$spec_var, .envs_a, "Site")
  expected <- 100 * sum(s$loads^2) / (sum(s$loads^2) + sum(.psi_a))
  expect_equal(v$vaf_total, expected, tolerance = 1e-10)
})

test_that("A: .fa_split_labels splits outer and inner level names", {
  lab <- .fa_split_labels(c("Site_Env01:Variety_Var01",
                            "Site_Env02:Variety_Var10"))
  expect_equal(lab$outer, c("Env01", "Env02"))
  expect_equal(lab$inner, c("Var01", "Var10"))
})

test_that("A: .fa_split_labels keeps underscores inside level names", {
  lab <- .fa_split_labels("Site_Env_01:Variety_Var_A_1")
  expect_equal(lab$outer, "Env_01")
  expect_equal(lab$inner, "Var_A_1")
})

test_that("A: .fa_split_labels handles Comp pseudo-levels", {
  lab <- .fa_split_labels("Site_Comp1:Variety_Var01")
  expect_equal(lab$outer, "Comp1")
})

# ===========================================================================
# SECTION B – Golden-value tests against real ASExtras4 output
# ===========================================================================

.golden_term <- function(fix_name, tm) {
  fx  <- .load_fix(fix_name)
  fas <- faSummary(fx$model)
  list(ref = fx$ref[[tm]], got = fas$gammas[[tm]], blup = fas$blups[[tm]])
}

for (.fx in c("fa_fix_k3", "fa_fix_k1", "fa_fix_vmide")) {
  local({
    fix_name <- .fx

    test_that(paste0("B: ", fix_name, " - Gmat and Cmat match ASExtras4"), {
      .skip_unless_fixture(fix_name)
      fx  <- .load_fix(fix_name)
      fas <- faSummary(fx$model)
      expect_setequal(fas$terms, names(fx$ref))
      for (tm in names(fx$ref)) {
        expect_equal(fas$gammas[[tm]]$Gmat, fx$ref[[tm]]$Gmat, tolerance = 1e-8)
        expect_equal(fas$gammas[[tm]]$Cmat, fx$ref[[tm]]$Cmat, tolerance = 1e-8)
      }
    })

    test_that(paste0("B: ", fix_name, " - loadings match, including sign"), {
      .skip_unless_fixture(fix_name)
      fx  <- .load_fix(fix_name)
      fas <- faSummary(fx$model)
      for (tm in names(fx$ref)) {
        if (is.null(fx$ref[[tm]]$loads)) next
        expect_equal(fas$gammas[[tm]]$loads, fx$ref[[tm]]$loads,
                     tolerance = 1e-8)
        expect_equal(fas$gammas[[tm]]$loads_cor, fx$ref[[tm]]$loads_cor,
                     tolerance = 1e-8)
      }
    })

    test_that(paste0("B: ", fix_name, " - specific variances match"), {
      .skip_unless_fixture(fix_name)
      fx  <- .load_fix(fix_name)
      fas <- faSummary(fx$model)
      for (tm in names(fx$ref)) {
        if (is.null(fx$ref[[tm]]$spec_var)) next
        expect_equal(unname(fas$gammas[[tm]]$spec_var),
                     unname(fx$ref[[tm]]$spec_var), tolerance = 1e-8)
      }
    })

    test_that(paste0("B: ", fix_name, " - VAF matches the ASExtras4 %vaf"), {
      .skip_unless_fixture(fix_name)
      fx  <- .load_fix(fix_name)
      fas <- faSummary(fx$model)
      for (tm in names(fx$ref)) {
        if (is.null(fx$ref[[tm]]$vaf_site)) next
        g <- fas$gammas[[tm]]
        mine <- 100 * as.matrix(g$vaf_env[, paste0("Factor", seq_len(g$k)),
                                          drop = FALSE])
        expect_equal(unname(mine),
                     unname(fx$ref[[tm]]$vaf_site[, seq_len(g$k), drop = FALSE]),
                     tolerance = 1e-8)
        expect_equal(g$vaf_total, fx$ref[[tm]]$vaf_total, tolerance = 1e-8)
      }
    })

    test_that(paste0("B: ", fix_name, " - BLUPs and scores match"), {
      .skip_unless_fixture(fix_name)
      fx  <- .load_fix(fix_name)
      fas <- faSummary(fx$model)
      for (tm in names(fx$ref)) {
        ref <- fx$ref[[tm]]$blups_in
        got <- fas$blups[[tm]]
        if (is.null(ref) || is.null(got)) next
        expect_equal(nrow(got$blups), nrow(ref))
        expect_equal(got$blups$blup, ref$blup, tolerance = 1e-8)
        expect_equal(got$blups$regblup, ref$regblup, tolerance = 1e-8)
        expect_equal(as.integer(got$blups$pres), as.integer(ref$pres))

        refs <- fx$ref[[tm]]$scores_in
        if (is.null(refs) || is.null(got$scores)) next
        expect_equal(got$scores$blup, refs$blup, tolerance = 1e-8)
        expect_equal(got$scores$blupr, refs$blupr, tolerance = 1e-8)
      }
    })
  })
}

test_that("B: the common variety effect is rotation invariant", {
  .skip_unless_fixture("fa_fix_k3")
  fx  <- .load_fix("fa_fix_k3")
  fas <- faSummary(fx$model)
  tm  <- fas$terms[1L]
  g   <- fas$gammas[[tm]]
  s   <- fas$blups[[tm]]$scores

  score_mat <- matrix(s$blupr, ncol = g$k)
  cve <- score_mat %*% t(g$loads)
  # regblup in the blups frame is exactly this, unrolled column-major
  expect_equal(as.vector(cve), fas$blups[[tm]]$blups$regblup,
               tolerance = 1e-8)
})

# ===========================================================================
# SECTION C – API and structure
# ===========================================================================

test_that("C: return value has the documented class and components", {
  .skip_unless_fixture("fa_fix_k3")
  fas <- faSummary(.load_fix("fa_fix_k3")$model)
  expect_s3_class(fas, "faSummary")
  expect_named(fas, c("gammas", "blups", "terms", "call"))
  g <- fas$gammas[[1L]]
  expect_true(all(c("Gmat", "Cmat", "loads", "loads_cor", "spec_var",
                    "vaf_env", "vaf_summary", "vaf_total", "k", "env",
                    "outer", "inner", "inner_fun") %in% names(g)))
})

test_that("C: blups and scores have the documented columns", {
  .skip_unless_fixture("fa_fix_k3")
  fas <- faSummary(.load_fix("fa_fix_k3")$model)
  b <- fas$blups[[1L]]
  expect_equal(names(b$blups), c("Site", "Variety", "blup", "regblup", "pres"))
  expect_equal(names(b$scores), c("Site", "Variety", "blup", "blupr"))
  expect_setequal(unique(b$scores$Site), paste0("Comp", 1:3))
})

test_that("C: term = selects a single FA term", {
  .skip_unless_fixture("fa_fix_vmide")
  m <- .load_fix("fa_fix_vmide")$model
  fas <- faSummary(m, term = "fa(Site, 2):ide(Variety)")
  expect_equal(fas$terms, "fa(Site, 2):ide(Variety)")
})

test_that("C: unknown term is rejected with the available terms listed", {
  .skip_unless_fixture("fa_fix_k3")
  m <- .load_fix("fa_fix_k3")$model
  expect_error(faSummary(m, term = "fa(Nope, 2):Variety"),
               "not found in the model")
})

test_that("C: blups = FALSE skips the BLUP extraction", {
  .skip_unless_fixture("fa_fix_k3")
  fas <- faSummary(.load_fix("fa_fix_k3")$model, blups = FALSE)
  expect_null(fas$blups)
  expect_false(is.null(fas$gammas[[1L]]$Gmat))
})

test_that("C: combine.ide appends the total structure", {
  .skip_unless_fixture("fa_fix_vmide")
  m <- .load_fix("fa_fix_vmide")$model
  fas <- faSummary(m, combine.ide = TRUE)
  expect_true("Site:Variety-total" %in% fas$terms)
  tot <- fas$gammas[["Site:Variety-total"]]
  expect_equal(tot$Gmat,
               fas$gammas[["fa(Site, 2):vm(Variety, Ginv)"]]$Gmat +
                 fas$gammas[["fa(Site, 2):ide(Variety)"]]$Gmat,
               tolerance = 1e-10)
  expect_equal(tot$inner_fun, "total")
})

test_that("C: combine.ide = FALSE suppresses the total structure", {
  .skip_unless_fixture("fa_fix_vmide")
  m <- .load_fix("fa_fix_vmide")$model
  fas <- faSummary(m, combine.ide = FALSE)
  expect_false("Site:Variety-total" %in% fas$terms)
  expect_length(fas$terms, 2L)
})

test_that("C: the combined term has no loadings or scores", {
  .skip_unless_fixture("fa_fix_vmide")
  fas <- faSummary(.load_fix("fa_fix_vmide")$model)
  expect_null(fas$gammas[["Site:Variety-total"]]$loads)
  expect_null(fas$blups[["Site:Variety-total"]]$scores)
})

test_that("C: k = 1 models are summarised without error", {
  .skip_unless_fixture("fa_fix_k1")
  fas <- faSummary(.load_fix("fa_fix_k1")$model)
  expect_equal(fas$gammas[[1L]]$k, 1L)
  expect_equal(ncol(fas$gammas[[1L]]$loads), 1L)
  expect_false(is.null(fas$blups[[1L]]$scores))
})

test_that("C: print method reports terms and variance accounted for", {
  .skip_unless_fixture("fa_fix_k3")
  fas <- faSummary(.load_fix("fa_fix_k3")$model)
  expect_output(print(fas), "faSummary")
  expect_output(print(fas), "k = 3")
  expect_output(print(fas), "variance accounted for")
})

# ===========================================================================
# SECTION D – Error handling
# ===========================================================================

test_that("D: a non-asreml object is rejected", {
  expect_error(faSummary(list(a = 1)), "must be a fitted 'asreml' object")
  expect_error(faSummary(data.frame(x = 1)), "must be a fitted 'asreml' object")
})

test_that("D: a model with no model frame and no RDS is rejected", {
  m <- structure(list(mf = NULL, RDS = NULL), class = "asreml")
  skip_if_not_installed("asreml")
  expect_error(faSummary(m), "missing its model frame")
})

test_that("D: .fa_find_terms errors when there is no random term", {
  mf <- data.frame(y = 1)
  expect_error(.fa_find_terms(mf), "no random terms")
})

test_that("D: .fa_find_terms errors when no fa() term is present", {
  mf <- data.frame(y = 1)
  tobj <- structure(list(), specials = list(fa = integer(0)))
  attr(mf, "model.terms") <- list(random = list(Terms.obj = tobj))
  expect_error(.fa_find_terms(mf), "No fa\\(\\) term found")
})
