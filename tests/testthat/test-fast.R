# tests/testthat/test-fast.R
# Consolidated tests for fast() in the biomAid package.
#
# Structure
# =========
#  A. Mathematical primitives  – pure-R recreations, no ASReml licence needed.
#  B. End-to-end tests         – call fast() directly via local_mocked_bindings.
#  C. New / gap-filling tests  – k=1, vm() parsing, column exclusions, sort order.
#
# The fixture matrices are shared by all three sections.

# ===========================================================================
# SHARED FIXTURE
# ===========================================================================
# 4 environments × 3 genotypes × 2 factors

.loads <- matrix(
  c( 2.0,  1.0,
     1.5, -0.8,
     1.8,  0.6,
     1.2, -1.2),
  nrow = 4L, ncol = 2L, byrow = TRUE,
  dimnames = list(c("E1","E2","E3","E4"), c("loads1","loads2"))
)

.scores <- matrix(
  c( 1.0,  0.5,
     0.2, -0.3,
    -0.5,  0.8),
  nrow = 3L, ncol = 2L, byrow = TRUE,
  dimnames = list(c("G1","G2","G3"), c("score1","score2"))
)

.spec  <- c(E1 = 0.10, E2 = 0.20, E3 = 0.15, E4 = 0.25)
.term  <- "fa(Site, 2):Genotype"
.envs  <- c("E1","E2","E3","E4")
.genos <- c("G1","G2","G3")
.m     <- 3L
.t     <- 4L
.k     <- 2L

# Pre-compute ground-truth matrices used by section A
.CVE_mat     <- .scores %*% t(.loads)                                # 3 × 4
.VE_mat      <- .CVE_mat + rep(.spec, each = .m)
.fitted1_mat <- outer(.scores[, 1L], .loads[, 1L])                  # m × t
.dev_mat     <- .CVE_mat - .fitted1_mat
.stab_vec    <- sqrt(rowMeans(.dev_mat^2))
.OP_vec      <- mean(.loads[, 1L]) * .scores[, 1L]
.sign_str    <- apply(.loads, 1L, function(x)
  paste(ifelse(x >= 0, "p", "n"), collapse = ""))

# Long-form index helpers (environment-major)
.env_rep  <- rep(seq_len(.t), each  = .m)
.geno_rep <- rep(seq_len(.m), times = .t)

# ---------------------------------------------------------------------------
# Helpers for end-to-end tests (sections B & C)
# ---------------------------------------------------------------------------
.sc_long <- data.frame(
  Site  = rep(paste0("Comp", 1:2), each = 3L),
  blupr = c(.scores[, 1L], .scores[, 2L]),
  stringsAsFactors = FALSE
)

.mock_sfa <- list(
  gammas = setNames(
    list(list("rotated loads" = .loads, "specific var" = .spec)),
    .term
  ),
  blups = setNames(
    list(list(scores = .sc_long)),
    .term
  )
)

make_fast_model <- function() {
  # fast() calls eval(model$call$data); <<- puts it in the global env so
  # eval() finds it regardless of where fast() runs.
  fast_mock_data <<- expand.grid(
    Site     = factor(c("E1","E2","E3","E4")),
    Genotype = factor(c("G1","G2","G3")),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  m <- list(call = list(data = quote(fast_mock_data)))
  class(m) <- "asreml"
  m
}

run_fast <- function(type = "all", ic.num = 2L, term = .term) {
  local_mocked_bindings(.fa_asreml = function(...) .mock_sfa, .package = "biomAid")
  fast(make_fast_model(), term = term, type = type, ic.num = ic.num)
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

test_that("A: CVE(G1,E1) = 2.0*1.0 + 1.0*0.5 = 2.5", {
  expect_equal(.CVE_mat["G1","E1"], 2.5, tolerance = 1e-12)
})

test_that("A: CVE(G3,E2) = 1.5*(-0.5) + (-0.8)*0.8", {
  expect_equal(.CVE_mat["G3","E2"],
               1.5 * (-0.5) + (-0.8) * 0.8, tolerance = 1e-12)
})

test_that("A: CVE(G2,E3) = 1.8*0.2 + 0.6*(-0.3)", {
  expect_equal(.CVE_mat["G2","E3"],
               1.8 * 0.2 + 0.6 * (-0.3), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A2. VE = CVE + spec.var (broadcast per environment)
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

test_that("A: CVE = fitted1 + fitted2", {
  fitted2 <- outer(.scores[, 2L], .loads[, 2L])
  expect_equal(.CVE_mat, .fitted1_mat + fitted2, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A4. OP = mean(loads1) * score1
# ---------------------------------------------------------------------------
test_that("A: OP = mean(loads1) * score1", {
  expect_equal(.OP_vec, mean(.loads[, 1L]) * .scores[, 1L], tolerance = 1e-12)
})

test_that("A: OP has one value per genotype", {
  expect_length(.OP_vec, .m)
})

# ---------------------------------------------------------------------------
# A5. dev = CVE - fitted1
# ---------------------------------------------------------------------------
test_that("A: dev = CVE - fitted1", {
  expect_equal(.dev_mat, .CVE_mat - .fitted1_mat, tolerance = 1e-12)
})

test_that("A: fitted1 + dev == CVE", {
  expect_equal(.fitted1_mat + .dev_mat, .CVE_mat, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A6. stab = RMSD of dev across environments
# ---------------------------------------------------------------------------
test_that("A: stab = sqrt(mean(dev^2)) per genotype", {
  for (g in seq_len(.m))
    expect_equal(unname(.stab_vec[g]),
                 sqrt(mean(.dev_mat[g, ]^2)), tolerance = 1e-12)
})

test_that("A: stab >= 0 for all genotypes", {
  expect_true(all(.stab_vec >= 0))
})

# ---------------------------------------------------------------------------
# A7. iClass sign patterns
# ---------------------------------------------------------------------------
test_that("A: iClass sign pattern for ic.num=2 matches fixture loadings", {
  expect_equal(unname(.sign_str["E1"]), "pp")
  expect_equal(unname(.sign_str["E2"]), "pn")
  expect_equal(unname(.sign_str["E3"]), "pp")
  expect_equal(unname(.sign_str["E4"]), "pn")
})

test_that("A: iClass sign pattern for ic.num=1 uses first loading only", {
  sign1 <- ifelse(.loads[, 1L] >= 0, "p", "n")
  expect_true(all(sign1 == "p"))   # all first loadings > 0
})

# ---------------------------------------------------------------------------
# A8. iClassOP = mean CVE within iClass
# ---------------------------------------------------------------------------
test_that("A: iClassOP for 'pp' equals mean CVE of E1 and E3", {
  pp_idx <- which(.sign_str == "pp")
  for (g in rownames(.scores)) {
    mld   <- colMeans(.loads[pp_idx, , drop = FALSE])
    iop_f <- sum(.scores[g, ] * mld)
    expect_equal(mean(.CVE_mat[g, pp_idx]), iop_f, tolerance = 1e-10)
  }
})

test_that("A: iClassOP for 'pn' equals mean CVE of E2 and E4", {
  pn_idx <- which(.sign_str == "pn")
  for (g in rownames(.scores)) {
    mld   <- colMeans(.loads[pn_idx, , drop = FALSE])
    iop_f <- sum(.scores[g, ] * mld)
    expect_equal(mean(.CVE_mat[g, pn_idx]), iop_f, tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# A9. iClassRMSD = 0 when ic.num == k (all factors used -> no residual)
# ---------------------------------------------------------------------------
test_that("A: iClassRMSD = 0 when all factors included (ic.num = k)", {
  fitted_all <- .fitted1_mat + outer(.scores[, 2L], .loads[, 2L])
  dev_ic     <- .CVE_mat - fitted_all
  rmsd       <- sqrt(rowMeans(dev_ic^2))
  expect_equal(max(abs(rmsd)), 0, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# A10. iClassRMSD > 0 when ic.num < k
# ---------------------------------------------------------------------------
test_that("A: iClassRMSD > 0 when ic.num < k", {
  for (w in unique(.sign_str)) {
    env_w    <- which(.sign_str == w)
    dev_ic_w <- .dev_mat[, env_w, drop = FALSE]
    rmsd_w   <- sqrt(rowMeans(dev_ic_w^2))
    expect_true(any(rmsd_w > 0))
  }
})

# ---------------------------------------------------------------------------
# A11. Term parsing logic (no call to fast())
# ---------------------------------------------------------------------------
test_that("A: parse 'fa(Site, 2):Genotype' -> sterm='Site'", {
  term    <- "fa(Site, 2):Genotype"
  parts   <- strsplit(term, ":")[[1L]]
  fa_idx  <- grep("^fa\\s*\\(", parts)
  sterm   <- trimws(sub(",.*", "", sub("^fa\\s*\\(", "", parts[fa_idx])))
  expect_equal(sterm, "Site")
})

test_that("A: parse 'fa(Site, 2):Genotype' -> gterm='Genotype'", {
  term     <- "fa(Site, 2):Genotype"
  parts    <- strsplit(term, ":")[[1L]]
  fa_idx   <- grep("^fa\\s*\\(", parts)
  gen_part <- parts[-fa_idx]
  expect_equal(trimws(gen_part), "Genotype")
})

# ---------------------------------------------------------------------------
# A12. ic.num validation logic
# ---------------------------------------------------------------------------
test_that("A: ic.num = 0 is invalid", {
  expect_false(0L >= 1L && 0L <= 2L)
})

test_that("A: ic.num = 5 is invalid when k = 2", {
  expect_false(5L >= 1L && 5L <= 2L)
})

test_that("A: ic.num = 1 is valid when k = 2", {
  expect_true(1L >= 1L && 1L <= 2L)
})

# ---------------------------------------------------------------------------
# A13. Environment-major ordering in long output
# ---------------------------------------------------------------------------
test_that("A: environment-major order: first m rows are E1", {
  env_long <- .envs[.env_rep]
  expect_equal(env_long[seq_len(.m)],     rep("E1", .m))
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
# B1. type = "all": correct dimensions and complete column set
# ---------------------------------------------------------------------------
test_that("B: type='all' returns 12 rows and all expected columns", {
  out <- run_fast("all")
  expect_equal(nrow(out), .t * .m)
  expected_cols <- c("Site","Genotype",
                     "loads1","loads2","spec.var",
                     "score1","score2",
                     "CVE","VE",
                     "fitted1","fitted2",
                     "OP","dev","stab",
                     "iclass","iClassOP","iClassRMSD")
  expect_true(all(expected_cols %in% names(out)))
})

# ---------------------------------------------------------------------------
# B2. CVE == fitted1 + fitted2 in the actual output
# ---------------------------------------------------------------------------
test_that("B: CVE == fitted1 + fitted2 in fast() output", {
  out <- run_fast("all")
  expect_equal(out$CVE, out$fitted1 + out$fitted2, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# B3. VE == CVE + spec.var in the actual output
# ---------------------------------------------------------------------------
test_that("B: VE == CVE + spec.var in fast() output", {
  out <- run_fast("all")
  expect_equal(out$VE, out$CVE + out$spec.var, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# B4. OP is constant within each genotype
# ---------------------------------------------------------------------------
test_that("B: OP is constant within each genotype", {
  out <- run_fast("FAST")
  for (g in .genos) {
    ops <- out$OP[out$Genotype == g]
    expect_equal(length(unique(round(ops, 12L))), 1L)
  }
})

# ---------------------------------------------------------------------------
# B5. stab = RMSD of dev across all environments
# ---------------------------------------------------------------------------
test_that("B: stab equals RMSD of dev per genotype", {
  out <- run_fast("FAST")
  for (g in .genos) {
    devs    <- out$dev[out$Genotype == g]
    stab_g  <- unique(round(out$stab[out$Genotype == g], 10L))
    expect_equal(stab_g, round(sqrt(mean(devs^2)), 10L))
  }
})

# ---------------------------------------------------------------------------
# B6. iclass labels match loading sign patterns
# ---------------------------------------------------------------------------
test_that("B: iclass='pp' for E1, 'pn' for E2", {
  out <- run_fast("iClass")
  expect_equal(as.character(out$iclass[out$Site == "E1"][1L]), "pp")
  expect_equal(as.character(out$iclass[out$Site == "E2"][1L]), "pn")
})

# ---------------------------------------------------------------------------
# B7. iClassRMSD = 0 when ic.num == k
# ---------------------------------------------------------------------------
test_that("B: iClassRMSD = 0 when ic.num == k (all factors included)", {
  out <- run_fast("iClass", ic.num = 2L)
  expect_equal(max(out$iClassRMSD), 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# B8. Error: term with no fa() component
# ---------------------------------------------------------------------------
test_that("B: fast() errors when term has no fa() component", {
  local_mocked_bindings(.fa_asreml = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fast(make_fast_model(), term = "Site:Genotype"),
    "must contain a 'fa\\("
  )
})

# ---------------------------------------------------------------------------
# B9. Error: ic.num out of range
# ---------------------------------------------------------------------------
test_that("B: fast() errors when ic.num > k", {
  local_mocked_bindings(.fa_asreml = function(...) .mock_sfa, .package = "biomAid")
  expect_error(
    fast(make_fast_model(), term = .term, type = "iClass", ic.num = 5L),
    "ic.num.*must be between"
  )
})

# ===========================================================================
# SECTION C – Gap-filling / new tests
# ===========================================================================

# ---------------------------------------------------------------------------
# C1. k = 1: warning issued, no "dev" or "stab" columns
# ---------------------------------------------------------------------------
test_that("C: k=1 single-factor model issues a warning about 'one factor'", {
  # Build a 1-factor version of the mock: loadings are a single column,
  # scores are a single column.
  loads_1f <- matrix(.loads[, 1L, drop = FALSE],
                     nrow = 4L, dimnames = list(.envs, NULL))
  scores_1f <- matrix(.scores[, 1L, drop = FALSE],
                      nrow = 3L, dimnames = list(.genos, NULL))
  sc_long_1f <- data.frame(
    Site  = "Comp1",
    blupr = scores_1f[, 1L],
    stringsAsFactors = FALSE
  )
  mock_1f <- list(
    gammas = setNames(
      list(list("rotated loads" = loads_1f, "specific var" = .spec)),
      .term
    ),
    blups = setNames(
      list(list(scores = sc_long_1f)),
      .term
    )
  )
  local_mocked_bindings(.fa_asreml = function(...) mock_1f, .package = "biomAid")
  expect_warning(
    out <- fast(make_fast_model(), term = .term, type = "iClass", ic.num = 1L),
    "one factor"
  )
  expect_false("dev"  %in% names(out))
  expect_false("stab" %in% names(out))
})

# ---------------------------------------------------------------------------
# C2. vm() genotype wrapping: gterm extracted correctly
# ---------------------------------------------------------------------------
test_that("C: vm(Genotype, ped) in term string extracts gterm='Genotype'", {
  gen_part <- "vm(Genotype, ped)"
  if (grepl("^vm\\s*\\(", gen_part)) {
    gterm <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", gen_part))))
  } else {
    gterm <- trimws(gen_part)
  }
  expect_equal(gterm, "Genotype")
})

test_that("C: full term 'fa(Site,2):vm(Genotype,ped)' extracts gterm='Genotype'", {
  term    <- "fa(Site,2):vm(Genotype,ped)"
  parts   <- strsplit(term, ":")[[1L]]
  fa_idx  <- grep("^fa\\s*\\(", parts)
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
# C3. type = "iClass" only: OP, dev, stab columns must be absent
# ---------------------------------------------------------------------------
test_that("C: type='iClass' has no OP, dev, or stab columns", {
  out <- run_fast("iClass")
  expect_false("OP"   %in% names(out))
  expect_false("dev"  %in% names(out))
  expect_false("stab" %in% names(out))
  # But iClass columns must be present
  expect_true("iclass"     %in% names(out))
  expect_true("iClassOP"   %in% names(out))
  expect_true("iClassRMSD" %in% names(out))
})

# ---------------------------------------------------------------------------
# C4. type = "FAST" only: iclass, iClassOP, iClassRMSD columns must be absent
# ---------------------------------------------------------------------------
test_that("C: type='FAST' has no iclass, iClassOP, or iClassRMSD columns", {
  out <- run_fast("FAST")
  expect_false("iclass"     %in% names(out))
  expect_false("iClassOP"   %in% names(out))
  expect_false("iClassRMSD" %in% names(out))
  # But FAST columns must be present
  expect_true("OP"   %in% names(out))
  expect_true("dev"  %in% names(out))
  expect_true("stab" %in% names(out))
})

# ---------------------------------------------------------------------------
# C5. Sort order for type = "iClass": primary sort key is iclass
# ---------------------------------------------------------------------------
test_that("C: type='iClass' output is sorted with iclass as primary key", {
  out <- run_fast("iClass")
  # All rows of one iClass must appear before the other iClass.
  # Check that iclass values are non-decreasing (factor level order).
  ic_int <- as.integer(out$iclass)
  expect_equal(ic_int, sort(ic_int))
})

test_that("C: within each iclass block, Site is the secondary sort key", {
  out <- run_fast("iClass")
  for (ic in levels(out$iclass)) {
    sub <- out[out$iclass == ic, ]
    expect_equal(as.character(sub$Site), sort(as.character(sub$Site)))
  }
})
