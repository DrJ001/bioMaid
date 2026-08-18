# tests/testthat/test-fast.R
# Consolidated tests for fastIC() in the biomAid package.
#
# Structure
# =========
#  A. Mathematical primitives  – pure-R recreations, no ASReml licence needed.
#  B. End-to-end tests         – call fastIC() directly via local_mocked_bindings.
#  C. Gap-filling tests        – k=1, vm()/ide() parsing, sort order.
#  D. ic.num boundary tests    – 1 <= ic.num < k enforced unconditionally.
#
# The fixture matrices are shared by all sections.
# NOTE: type argument removed — fastIC() always computes both global_* and
# iClass* columns. ic.num validation is now unconditional.

# ===========================================================================
# SHARED FIXTURE
# ===========================================================================
# 4 environments x 3 genotypes x 3 factors
# ic.num = 2 (default) is valid: 2 < k = 3.
# Factor 3 is reserved for iClassRMSD.

.loads <- matrix(
  c( 2.0,  1.0,  0.4,
     1.5, -0.8,  0.6,
     1.8,  0.6, -0.5,
     1.2, -1.2,  0.3),
  nrow = 4L, ncol = 3L, byrow = TRUE,
  dimnames = list(c("E1","E2","E3","E4"), c("loads1","loads2","loads3"))
)

.scores <- matrix(
  c( 1.0,  0.5,  0.3,
     0.2, -0.3,  0.4,
    -0.5,  0.8, -0.2),
  nrow = 3L, ncol = 3L, byrow = TRUE,
  dimnames = list(c("G1","G2","G3"), c("score1","score2","score3"))
)

.spec  <- c(E1 = 0.10, E2 = 0.20, E3 = 0.15, E4 = 0.25)
.term  <- "fa(Site, 3):Genotype"
.envs  <- c("E1","E2","E3","E4")
.genos <- c("G1","G2","G3")
.m     <- 3L
.t     <- 4L
.k     <- 3L

# Pre-compute ground-truth matrices
.CVE_mat     <- .scores %*% t(.loads)
.VE_mat      <- .CVE_mat + rep(.spec, each = .m)
.fitted1_mat <- outer(.scores[, 1L], .loads[, 1L])
.dev_mat     <- .CVE_mat - .fitted1_mat
.stab_vec    <- sqrt(rowMeans(.dev_mat^2))
.OP_vec      <- mean(.loads[, 1L]) * .scores[, 1L]

# Sign pattern for ic.num = 2 (first 2 factors)
.sign_str <- apply(.loads[, 1:2, drop = FALSE], 1L, function(x)
  paste(ifelse(x >= 0, "p", "n"), collapse = ""))

# Long-form index helpers (environment-major)
.env_rep  <- rep(seq_len(.t), each  = .m)
.geno_rep <- rep(seq_len(.m), times = .t)

# ---------------------------------------------------------------------------
# Helpers for end-to-end tests
# ---------------------------------------------------------------------------
.sc_long <- data.frame(
  Site  = rep(paste0("Comp", 1:3), each = 3L),
  blupr = c(.scores[, 1L], .scores[, 2L], .scores[, 3L]),
  stringsAsFactors = FALSE
)

# Build a faSummary()-shaped mock.  VAF is produced by the real internal
# helper so the mock stays consistent with the production calculation.
.mk_mock_sfa <- function(loads, spec, sc_long, term, envs = .envs) {
  vaf <- .fa_vaf(loads, spec, envs, "Site")
  structure(
    list(
      gammas = setNames(
        list(list(
          loads       = loads,
          spec_var    = spec,
          vaf_env     = vaf$vaf_env,
          vaf_summary = vaf$vaf_summary,
          vaf_total   = vaf$vaf_total,
          k           = ncol(loads),
          env         = envs
        )),
        term
      ),
      blups = setNames(list(list(scores = sc_long)), term),
      terms = term
    ),
    class = "faSummary"
  )
}

.mock_sfa <- .mk_mock_sfa(.loads, .spec, .sc_long, .term)

make_fastIC_model <- function() {
  fastIC_mock_data <<- expand.grid(
    Site     = factor(c("E1","E2","E3","E4")),
    Genotype = factor(c("G1","G2","G3")),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  m <- list(call = list(data = quote(fastIC_mock_data)))
  class(m) <- "asreml"
  m
}

run_fastIC <- function(ic.num = 2L, term = .term) {
  local_mocked_bindings(faSummary = function(...) .mock_sfa, .package = "biomAid")
  fastIC(make_fastIC_model(), term = term, ic.num = ic.num)
}

# ===========================================================================
# SECTION A – Mathematical primitives (pure R, no ASReml)
# ===========================================================================

# ---------------------------------------------------------------------------
# A1. CVE = score_mat %*% t(loads_mat)
# ---------------------------------------------------------------------------
test_that("A: CVE matrix has correct dimensions", {
  expect_equal(dim(.CVE_mat), c(.m, .t))
})

test_that("A: CVE(G1,E1) = 2.0*1.0 + 1.0*0.5 + 0.4*0.3", {
  expect_equal(.CVE_mat["G1","E1"],
               2.0 * 1.0 + 1.0 * 0.5 + 0.4 * 0.3, tolerance = 1e-12)
})

test_that("A: CVE(G3,E2) = 1.5*(-0.5) + (-0.8)*0.8 + 0.6*(-0.2)", {
  expect_equal(.CVE_mat["G3","E2"],
               1.5 * (-0.5) + (-0.8) * 0.8 + 0.6 * (-0.2), tolerance = 1e-12)
})

test_that("A: CVE(G2,E3) = 1.8*0.2 + 0.6*(-0.3) + (-0.5)*0.4", {
  expect_equal(.CVE_mat["G2","E3"],
               1.8 * 0.2 + 0.6 * (-0.3) + (-0.5) * 0.4, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A2. VE = CVE + spec.var
# ---------------------------------------------------------------------------
test_that("A: VE = CVE + spec.var for each environment column", {
  for (j in seq_len(.t))
    expect_equal(.VE_mat[, j], .CVE_mat[, j] + .spec[j], tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A3. Fitted factor outer products
# ---------------------------------------------------------------------------
test_that("A: fitted1 = outer(score1, loads1)", {
  expect_equal(.fitted1_mat, outer(.scores[, 1L], .loads[, 1L]), tolerance = 1e-12)
})

test_that("A: CVE = fitted1 + fitted2 + fitted3", {
  fitted2 <- outer(.scores[, 2L], .loads[, 2L])
  fitted3 <- outer(.scores[, 3L], .loads[, 3L])
  expect_equal(.CVE_mat, .fitted1_mat + fitted2 + fitted3, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A4. global_op = mean(loads1) * score1
# ---------------------------------------------------------------------------
test_that("A: global_op = mean(loads1) * score1", {
  expect_equal(.OP_vec, mean(.loads[, 1L]) * .scores[, 1L], tolerance = 1e-12)
})

test_that("A: global_op has one value per genotype", {
  expect_length(.OP_vec, .m)
})

# ---------------------------------------------------------------------------
# A5. global_dev = CVE - fitted1
# ---------------------------------------------------------------------------
test_that("A: global_dev = CVE - fitted1", {
  expect_equal(.dev_mat, .CVE_mat - .fitted1_mat, tolerance = 1e-12)
})

test_that("A: fitted1 + global_dev == CVE", {
  expect_equal(.fitted1_mat + .dev_mat, .CVE_mat, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A6. global_stab = RMSD of global_dev across environments
# ---------------------------------------------------------------------------
test_that("A: global_stab = sqrt(mean(global_dev^2)) per genotype", {
  for (g in seq_len(.m))
    expect_equal(unname(.stab_vec[g]),
                 sqrt(mean(.dev_mat[g, ]^2)), tolerance = 1e-12)
})

test_that("A: global_stab >= 0 for all genotypes", {
  expect_true(all(.stab_vec >= 0))
})

# ---------------------------------------------------------------------------
# A7. iClass sign patterns (ic.num = 2: first 2 factors)
# E1: loads1=2.0(p), loads2=1.0(p)  -> "pp"
# E2: loads1=1.5(p), loads2=-0.8(n) -> "pn"
# E3: loads1=1.8(p), loads2=0.6(p)  -> "pp"
# E4: loads1=1.2(p), loads2=-1.2(n) -> "pn"
# ---------------------------------------------------------------------------
test_that("A: iClass sign pattern for ic.num=2 matches first 2 factor loadings", {
  expect_equal(unname(.sign_str["E1"]), "pp")
  expect_equal(unname(.sign_str["E2"]), "pn")
  expect_equal(unname(.sign_str["E3"]), "pp")
  expect_equal(unname(.sign_str["E4"]), "pn")
})

test_that("A: iClass sign pattern for ic.num=1 uses first loading only", {
  sign1 <- ifelse(.loads[, 1L] >= 0, "p", "n")
  expect_true(all(sign1 == "p"))
})

# ---------------------------------------------------------------------------
# A8. iClassOP formula
# ---------------------------------------------------------------------------
test_that("A: iClassOP formula: score[1:ic.num] dot mean_loads[1:ic.num] within class", {
  ic.num <- 2L
  for (w in unique(.sign_str)) {
    env_w <- which(.sign_str == w)
    mld_w <- colMeans(.loads[env_w, seq_len(ic.num), drop = FALSE])
    for (g in rownames(.scores)) {
      expected <- sum(.scores[g, seq_len(ic.num)] * mld_w)
      mat_form <- unname((.scores[g, seq_len(ic.num), drop = FALSE] %*% mld_w)[1, 1])
      expect_equal(expected, mat_form, tolerance = 1e-12)
    }
  }
})

test_that("A: iClassOP != mean(CVE) within class when ic.num < k", {
  ic.num <- 2L
  for (w in unique(.sign_str)) {
    env_w      <- which(.sign_str == w)
    mld_w      <- colMeans(.loads[env_w, seq_len(ic.num), drop = FALSE])
    iop_g      <- as.vector(.scores[, seq_len(ic.num), drop = FALSE] %*% mld_w)
    mean_cve_g <- rowMeans(.CVE_mat[, env_w, drop = FALSE])
    expect_false(isTRUE(all.equal(iop_g, mean_cve_g, tolerance = 1e-6)))
  }
})

# ---------------------------------------------------------------------------
# A9. iClassRMSD > 0 when ic.num = k-1
# ---------------------------------------------------------------------------
test_that("A: iClassRMSD > 0 when ic.num = k-1 (one factor reserved)", {
  fitted_ic2 <- outer(.scores[, 1L], .loads[, 1L]) +
                outer(.scores[, 2L], .loads[, 2L])
  dev_ic <- .CVE_mat - fitted_ic2
  rmsd   <- sqrt(rowMeans(dev_ic^2))
  expect_true(any(rmsd > 0))
})

# ---------------------------------------------------------------------------
# A10. iClassRMSD > 0 within each iClass
# ---------------------------------------------------------------------------
test_that("A: iClassRMSD > 0 within each iClass when ic.num < k", {
  fitted_ic2 <- outer(.scores[, 1L], .loads[, 1L]) +
                outer(.scores[, 2L], .loads[, 2L])
  for (w in unique(.sign_str)) {
    env_w    <- which(.sign_str == w)
    dev_ic_w <- (.CVE_mat - fitted_ic2)[, env_w, drop = FALSE]
    rmsd_w   <- sqrt(rowMeans(dev_ic_w^2))
    expect_true(any(rmsd_w > 0))
  }
})

# ---------------------------------------------------------------------------
# A11. Term parsing
# ---------------------------------------------------------------------------
test_that("A: parse 'fa(Site, 3):Genotype' -> sterm='Site'", {
  term   <- "fa(Site, 3):Genotype"
  parts  <- strsplit(term, ":")[[1L]]
  fa_idx <- grep("^fa\\s*\\(", parts)
  sterm  <- trimws(sub(",.*", "", sub("^fa\\s*\\(", "", parts[fa_idx])))
  expect_equal(sterm, "Site")
})

test_that("A: parse 'fa(Site, 3):Genotype' -> gterm='Genotype'", {
  term     <- "fa(Site, 3):Genotype"
  parts    <- strsplit(term, ":")[[1L]]
  fa_idx   <- grep("^fa\\s*\\(", parts)
  gen_part <- parts[-fa_idx]
  expect_equal(trimws(gen_part), "Genotype")
})

# ---------------------------------------------------------------------------
# A12. ic.num validation (1 <= ic.num < k, enforced unconditionally)
# ---------------------------------------------------------------------------
test_that("A: ic.num = 0 is invalid", {
  k <- 3L
  expect_false(0L >= 1L && 0L < k)
})

test_that("A: ic.num = k is invalid (must be strictly less than k)", {
  k <- 3L
  expect_false(3L >= 1L && 3L < k)
})

test_that("A: ic.num = k-1 is valid (upper boundary)", {
  k <- 3L
  expect_true(2L >= 1L && 2L < k)
})

test_that("A: ic.num = 1 is valid when k = 3", {
  k <- 3L
  expect_true(1L >= 1L && 1L < k)
})

# ---------------------------------------------------------------------------
# A13. Environment-major ordering
# ---------------------------------------------------------------------------
test_that("A: environment-major order: first m rows are E1", {
  env_long <- .envs[.env_rep]
  expect_equal(env_long[seq_len(.m)],         rep("E1", .m))
  expect_equal(env_long[(.m + 1L):(2L * .m)], rep("E2", .m))
})

# ---------------------------------------------------------------------------
# A14. spec.var broadcast
# ---------------------------------------------------------------------------
test_that("A: spec.var repeated correctly in long format", {
  sv_long <- .spec[.env_rep]
  for (e in .envs) {
    idx <- which(.envs[.env_rep] == e)
    expect_equal(length(unique(sv_long[idx])), 1L)
    expect_equal(unname(unique(sv_long[idx])), unname(.spec[e]))
  }
})

# ---------------------------------------------------------------------------
# A15. Matrix dimensions
# ---------------------------------------------------------------------------
test_that("A: CVE_mat is m x t", {
  expect_equal(dim(.CVE_mat), c(.m, .t))
})

test_that("A: dev_mat matches CVE_mat dimensions", {
  expect_equal(dim(.dev_mat), dim(.CVE_mat))
})

test_that("A: stab_vec has length m", {
  expect_length(.stab_vec, .m)
})

# ===========================================================================
# SECTION B – End-to-end tests via local_mocked_bindings
# ===========================================================================

# ---------------------------------------------------------------------------
# B1. Correct dimensions and complete column set (always all columns)
# ---------------------------------------------------------------------------
test_that("B: fastIC() returns 12 rows and all expected columns", {
  out <- run_fastIC()
  expect_equal(nrow(out), .t * .m)
  expected_cols <- c("Site", "Genotype",
                     "loads1", "loads2", "loads3", "spec.var",
                     "score1", "score2", "score3",
                     "CVE",
                     "fitted1", "fitted2", "fitted3",
                     "global_op", "global_dev", "global_stab",
                     "iclass", "iClassOP", "iClassRMSD")
  expect_true(all(expected_cols %in% names(out)))
})

# ---------------------------------------------------------------------------
# B2. CVE == fitted1 + fitted2 + fitted3
# ---------------------------------------------------------------------------
test_that("B: CVE == fitted1 + fitted2 + fitted3", {
  out <- run_fastIC()
  expect_equal(out$CVE, out$fitted1 + out$fitted2 + out$fitted3,
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# B3. VAF attributes are attached to the result
# ---------------------------------------------------------------------------
test_that("B: fastIC() result has vaf_env and vaf_summary attributes", {
  out <- run_fastIC()
  expect_false(is.null(attr(out, "vaf_env")))
  expect_false(is.null(attr(out, "vaf_summary")))
})

test_that("B: vaf_env has one row per environment with proportions summing to 1", {
  out     <- run_fastIC()
  vaf_env <- attr(out, "vaf_env")
  expect_equal(nrow(vaf_env), length(.envs))
  fac_cols <- paste0("Factor", 1:.k)
  row_sums <- rowSums(vaf_env[, c(fac_cols, "Specific"), drop = FALSE])
  expect_true(all(abs(row_sums - 1.0) < 1e-10))
})

test_that("B: vaf_summary has k+1 rows and pct_var sums to 1", {
  out  <- run_fastIC()
  summ <- attr(out, "vaf_summary")
  expect_equal(nrow(summ), .k + 1L)
  expect_equal(sum(summ$pct_var), 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# B4. global_op is constant within each genotype
# ---------------------------------------------------------------------------
test_that("B: global_op is constant within each genotype", {
  out <- run_fastIC()
  for (g in .genos) {
    ops <- out$global_op[out$Genotype == g]
    expect_equal(length(unique(round(ops, 12L))), 1L)
  }
})

# ---------------------------------------------------------------------------
# B5. global_stab = RMSD of global_dev per genotype
# ---------------------------------------------------------------------------
test_that("B: global_stab equals RMSD of global_dev per genotype", {
  out <- run_fastIC()
  for (g in .genos) {
    devs   <- out$global_dev[out$Genotype == g]
    stab_g <- unique(round(out$global_stab[out$Genotype == g], 10L))
    expect_equal(stab_g, round(sqrt(mean(devs^2)), 10L))
  }
})

# ---------------------------------------------------------------------------
# B6. global_op == mean(loads1) * score1 in output
# ---------------------------------------------------------------------------
test_that("B: global_op equals mean(loads1) * score1", {
  out    <- run_fastIC()
  exp_op <- mean(.loads[, 1L]) * out$score1
  expect_equal(out$global_op, exp_op, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# B7. global_dev == CVE - fitted1 in output
# ---------------------------------------------------------------------------
test_that("B: global_dev == CVE - fitted1", {
  out <- run_fastIC()
  expect_equal(out$global_dev, out$CVE - out$fitted1, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# B8. iclass labels match sign pattern of first ic.num=2 loadings
# ---------------------------------------------------------------------------
test_that("B: iclass='pp' for E1, 'pn' for E2 (based on factors 1 and 2)", {
  out <- run_fastIC(ic.num = 2L)
  expect_equal(as.character(out$iclass[out$Site == "E1"][1L]), "pp")
  expect_equal(as.character(out$iclass[out$Site == "E2"][1L]), "pn")
})

# ---------------------------------------------------------------------------
# B9. iClassRMSD > 0 when ic.num = k-1
# ---------------------------------------------------------------------------
test_that("B: iClassRMSD > 0 when ic.num = k-1 (kth factor reserved)", {
  out <- run_fastIC(ic.num = 2L)
  expect_true(any(out$iClassRMSD > 0))
})

# ---------------------------------------------------------------------------
# B10. Error: term with no fa() component
# ---------------------------------------------------------------------------
test_that("B: fastIC() errors when term has no fa() component", {
  local_mocked_bindings(faSummary = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fastIC(make_fastIC_model(), term = "Site:Genotype"),
    "must contain a 'fa\\("
  )
})

# ---------------------------------------------------------------------------
# B11. Error: ic.num >= k
# ---------------------------------------------------------------------------
test_that("B: fastIC() errors when ic.num = k", {
  local_mocked_bindings(faSummary = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fastIC(make_fastIC_model(), term = .term, ic.num = 3L),
    "ic.num.*must be between"
  )
})

test_that("B: fastIC() errors when ic.num > k", {
  local_mocked_bindings(faSummary = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fastIC(make_fastIC_model(), term = .term, ic.num = 5L),
    "ic.num.*must be between"
  )
})

# ---------------------------------------------------------------------------
# B12. When ic.num=1, iClassOP converges with global_op (single class)
#      All first-factor loadings in fixture are positive -> one class "p"
# ---------------------------------------------------------------------------
test_that("B: ic.num=1 produces one iClass and iClassOP == global_op", {
  out <- run_fastIC(ic.num = 1L)
  # Only one iClass level
  expect_equal(length(levels(out$iclass)), 1L)
  expect_equal(levels(out$iclass), "p")
  # iClassOP should equal global_op (within floating point)
  expect_equal(out$iClassOP, out$global_op, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# B13. ic.num=1: iClassRMSD == global_stab (single class, all environments)
# ---------------------------------------------------------------------------
test_that("B: ic.num=1 iClassRMSD equals global_stab", {
  out <- run_fastIC(ic.num = 1L)
  expect_equal(out$iClassRMSD, out$global_stab, tolerance = 1e-10)
})

# ===========================================================================
# SECTION C – Gap-filling / new tests
# ===========================================================================

# ---------------------------------------------------------------------------
# C1. k=1: warning issued, then ic.num validation errors (always runs now)
# ---------------------------------------------------------------------------
test_that("C: k=1 issues warning then errors on ic.num (validation always runs)", {
  loads_1f   <- matrix(.loads[, 1L, drop = FALSE],
                       nrow = 4L, dimnames = list(.envs, NULL))
  sc_long_1f <- data.frame(Site  = "Comp1",
                            blupr = .scores[, 1L],
                            stringsAsFactors = FALSE)
  mock_1f <- .mk_mock_sfa(loads_1f, .spec, sc_long_1f, .term)
  local_mocked_bindings(faSummary = function(...) mock_1f, .package = "biomAid")
  # warning fires (k==1), then stop fires (ic.num=1 >= k=1)
  expect_error(
    withCallingHandlers(
      fastIC(make_fastIC_model(), term = .term, ic.num = 1L),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    "ic.num.*must be between"
  )
})

# ---------------------------------------------------------------------------
# C2. vm() genotype wrapping: gterm extracted correctly
# ---------------------------------------------------------------------------
test_that("C: vm(Genotype, ped) extracts gterm='Genotype'", {
  gen_part <- "vm(Genotype, ped)"
  if (grepl("^vm\\s*\\(", gen_part)) {
    gterm <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", gen_part))))
  } else {
    gterm <- trimws(gen_part)
  }
  expect_equal(gterm, "Genotype")
})

test_that("C: full term 'fa(Site,3):vm(Genotype,ped)' extracts gterm='Genotype'", {
  term     <- "fa(Site,3):vm(Genotype,ped)"
  parts    <- strsplit(term, ":")[[1L]]
  fa_idx   <- grep("^fa\\s*\\(", parts)
  gen_part <- parts[-fa_idx]
  if (any(grepl("^vm\\s*\\(", gen_part))) {
    vm_part <- gen_part[grep("^vm\\s*\\(", gen_part)]
    gterm   <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", vm_part))))
  } else {
    gterm <- trimws(gen_part)
  }
  expect_equal(gterm, "Genotype")
})

# ---------------------------------------------------------------------------
# C3. ide() genotype wrapping: gterm extracted correctly
# ---------------------------------------------------------------------------
test_that("C: ide(Genotype) extracts gterm='Genotype'", {
  gen_part <- "ide(Genotype)"
  if (grepl("^ide\\s*\\(", gen_part)) {
    gterm <- trimws(sub("[,)].*", "", sub("^ide\\s*\\(", "", gen_part)))
  } else {
    gterm <- trimws(gen_part)
  }
  expect_equal(gterm, "Genotype")
})

test_that("C: full term 'fa(Site,3):ide(Genotype)' extracts gterm='Genotype'", {
  term     <- "fa(Site,3):ide(Genotype)"
  parts    <- strsplit(term, ":")[[1L]]
  fa_idx   <- grep("^fa\\s*\\(", parts)
  gen_part <- parts[-fa_idx]
  if (any(grepl("^vm\\s*\\(", gen_part))) {
    vm_part <- gen_part[grep("^vm\\s*\\(", gen_part)]
    gterm   <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", vm_part))))
  } else if (any(grepl("^ide\\s*\\(", gen_part))) {
    ide_part <- gen_part[grep("^ide\\s*\\(", gen_part)]
    gterm    <- trimws(sub("[,)].*", "", sub("^ide\\s*\\(", "", ide_part)))
  } else {
    gterm <- trimws(gen_part)
  }
  expect_equal(gterm, "Genotype")
})

test_that("C: fastIC() end-to-end works with ide() wrapper in term", {
  ide_term <- "fa(Site, 3):ide(Genotype)"
  mock_ide <- .mk_mock_sfa(.loads, .spec, .sc_long, ide_term)
  local_mocked_bindings(faSummary = function(...) mock_ide, .package = "biomAid")
  out <- fastIC(make_fastIC_model(), term = ide_term, ic.num = 2L)
  expect_true("Genotype"   %in% names(out))
  expect_true("global_op"  %in% names(out))
  expect_true("iclass"     %in% names(out))
  expect_equal(nrow(out), .t * .m)
})

# ---------------------------------------------------------------------------
# C4. Output always contains both global_* and iClass* columns
# ---------------------------------------------------------------------------
test_that("C: output always has global_op, global_stab, iclass, iClassOP, iClassRMSD", {
  out <- run_fastIC()
  for (col in c("global_op", "global_dev", "global_stab",
                "iclass", "iClassOP", "iClassRMSD"))
    expect_true(col %in% names(out),
                label = paste("column", col, "present"))
})

# ---------------------------------------------------------------------------
# C5. Sort order: primary = iclass, secondary = Site, tertiary = Genotype
# ---------------------------------------------------------------------------
test_that("C: output is sorted with iclass as primary key", {
  out    <- run_fastIC()
  ic_int <- as.integer(out$iclass)
  expect_equal(ic_int, sort(ic_int))
})

test_that("C: within each iclass block, Site is the secondary sort key", {
  out <- run_fastIC()
  for (ic in levels(out$iclass)) {
    sub <- out[out$iclass == ic, ]
    expect_equal(as.character(sub$Site), sort(as.character(sub$Site)))
  }
})

# ===========================================================================
# SECTION D – ic.num boundary validation (now always enforced)
# ===========================================================================

# ---------------------------------------------------------------------------
# D1. ic.num = k-1 is the valid upper boundary
# ---------------------------------------------------------------------------
test_that("D: ic.num = k-1 is accepted (valid upper boundary)", {
  expect_no_error(run_fastIC(ic.num = 2L))
})

test_that("D: ic.num = 1 is accepted (valid lower boundary)", {
  expect_no_error(run_fastIC(ic.num = 1L))
})

# ---------------------------------------------------------------------------
# D2. ic.num = k is invalid
# ---------------------------------------------------------------------------
test_that("D: ic.num = k errors with informative message", {
  local_mocked_bindings(faSummary = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fastIC(make_fastIC_model(), term = .term, ic.num = 3L),
    "kth factor must remain"
  )
})

# ---------------------------------------------------------------------------
# D3. ic.num = 0 is invalid
# ---------------------------------------------------------------------------
test_that("D: ic.num = 0 errors", {
  local_mocked_bindings(faSummary = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fastIC(make_fastIC_model(), term = .term, ic.num = 0L),
    "ic.num.*must be between"
  )
})

# ---------------------------------------------------------------------------
# D4. iClassRMSD values are non-negative
# ---------------------------------------------------------------------------
test_that("D: all iClassRMSD values are >= 0 for ic.num=1", {
  out <- run_fastIC(ic.num = 1L)
  expect_true(all(out$iClassRMSD >= 0))
})

test_that("D: all iClassRMSD values are >= 0 at upper boundary ic.num=k-1", {
  out <- run_fastIC(ic.num = 2L)
  expect_true(all(out$iClassRMSD >= 0))
})

# ---------------------------------------------------------------------------
# D5. global_stab values are non-negative
# ---------------------------------------------------------------------------
test_that("D: all global_stab values are >= 0", {
  out <- run_fastIC()
  expect_true(all(out$global_stab >= 0))
})
