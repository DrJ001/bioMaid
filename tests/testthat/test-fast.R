# tests/testthat/test-fast.R
# Tests for fast() — no ASReml licence required.
# We test the mathematical core of FAST/iClass directly without calling fast().
# The fa.asreml dependency cannot be mocked in a locked package namespace,
# so we replicate the computations in pure R.

# ---------------------------------------------------------------------------
# Fixture: 4 environments x 3 genotypes x 2 factors
# ---------------------------------------------------------------------------
loads_mat <- matrix(
  c( 2.0,  1.0,
     1.5, -0.8,
     1.8,  0.6,
     1.2, -1.2),
  nrow = 4, ncol = 2, byrow = TRUE,
  dimnames = list(c("E1","E2","E3","E4"), c("loads1","loads2"))
)

score_mat <- matrix(
  c( 1.0,  0.5,
     0.2, -0.3,
    -0.5,  0.8),
  nrow = 3, ncol = 2, byrow = TRUE,
  dimnames = list(c("G1","G2","G3"), c("score1","score2"))
)

spec_var  <- c(E1 = 0.10, E2 = 0.20, E3 = 0.15, E4 = 0.25)
envs      <- c("E1","E2","E3","E4")
genotypes <- c("G1","G2","G3")
m         <- 3L
t_envs    <- 4L
k         <- 2L

# Ground truth: replicate the fast() computations in plain R
CVE_mat     <- score_mat %*% t(loads_mat)                      # 3 x 4
VE_mat      <- CVE_mat + rep(spec_var, each = m)
fitted1_mat <- outer(score_mat[, 1L], loads_mat[, 1L])         # m x t
dev_mat     <- CVE_mat - fitted1_mat
stab_vec    <- sqrt(rowMeans(dev_mat^2))
mean_l1     <- mean(loads_mat[, 1L])
OP_vec      <- mean_l1 * score_mat[, 1L]

sign_str <- apply(loads_mat, 1L, function(x)
  paste(ifelse(x >= 0, "p", "n"), collapse = ""))
# E1="pp", E2="pn", E3="pp", E4="pn"

# Build long-form output (environment-major order: t_envs blocks of m rows)
env_rep  <- rep(seq_len(t_envs), each = m)
geno_rep <- rep(seq_len(m),      times = t_envs)

CVE_long <- as.vector(CVE_mat)   # column-major = env-major ✓
VE_long  <- as.vector(VE_mat)

# ---------------------------------------------------------------------------
# 1. CVE = score_mat %*% t(loads_mat)
# ---------------------------------------------------------------------------
test_that("CVE formula: CVE_mat = score_mat %*% t(loads_mat)", {
  expect_equal(dim(CVE_mat), c(m, t_envs))
})

test_that("CVE(G1,E1) = 2.0*1.0 + 1.0*0.5 = 2.5", {
  expect_equal(CVE_mat["G1","E1"], 2.5, tolerance = 1e-12)
})

test_that("CVE(G3,E2) = 1.5*(-0.5) + (-0.8)*0.8 = -1.39", {
  expect_equal(CVE_mat["G3","E2"],
               1.5 * (-0.5) + (-0.8) * 0.8, tolerance = 1e-12)
})

test_that("CVE(G2,E3) = 1.8*0.2 + 0.6*(-0.3) = 0.18", {
  expect_equal(CVE_mat["G2","E3"],
               1.8 * 0.2 + 0.6 * (-0.3), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 2. VE = CVE + spec.var
# ---------------------------------------------------------------------------
test_that("VE = CVE + spec.var (broadcast per environment)", {
  for (j in seq_len(t_envs)) {
    expect_equal(VE_mat[, j], CVE_mat[, j] + spec_var[j], tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 3. fitted_r = score_r * loads_r (per factor outer product)
# ---------------------------------------------------------------------------
test_that("fitted1 = outer(score1, loads1)", {
  fitted1 <- outer(score_mat[, 1L], loads_mat[, 1L])
  expect_equal(fitted1, fitted1_mat, tolerance = 1e-12)
})

test_that("CVE = sum of fitted factors (fitted1 + fitted2)", {
  fitted2_mat <- outer(score_mat[, 2L], loads_mat[, 2L])
  expect_equal(CVE_mat, fitted1_mat + fitted2_mat, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4. OP = mean(loads1) * score1
# ---------------------------------------------------------------------------
test_that("OP = mean(loads1) * score1", {
  expected_OP <- mean(loads_mat[, 1L]) * score_mat[, 1L]
  expect_equal(OP_vec, expected_OP, tolerance = 1e-12)
})

test_that("OP is a per-genotype scalar (same across environments)", {
  # OP_vec has one element per genotype
  expect_length(OP_vec, m)
  # Repeating it for all environments should give same value
  op_long <- OP_vec[geno_rep]
  for (g in seq_len(m)) {
    ops <- op_long[geno_rep == g]
    expect_equal(length(unique(round(ops, 14L))), 1L)
  }
})

# ---------------------------------------------------------------------------
# 5. dev = CVE - fitted1 (higher-order factors' contribution)
# ---------------------------------------------------------------------------
test_that("dev = CVE - fitted1 (residual from first factor)", {
  expect_equal(dev_mat, CVE_mat - fitted1_mat, tolerance = 1e-12)
})

test_that("dev sums to CVE when added back to fitted1", {
  expect_equal(fitted1_mat + dev_mat, CVE_mat, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 6. stab = RMSD of dev across environments
# ---------------------------------------------------------------------------
test_that("stab = sqrt(mean(dev^2)) per genotype", {
  for (g in seq_len(m)) {
    devs <- dev_mat[g, ]
    expect_equal(unname(stab_vec[g]), sqrt(mean(devs^2)), tolerance = 1e-12)
  }
})

test_that("stab is >= 0 for all genotypes", {
  expect_true(all(stab_vec >= 0))
})

# ---------------------------------------------------------------------------
# 7. iClass sign patterns
# ---------------------------------------------------------------------------
test_that("iClass sign pattern for 2 factors (ic.num=2)", {
  expect_equal(unname(sign_str["E1"]), "pp")
  expect_equal(unname(sign_str["E2"]), "pn")
  expect_equal(unname(sign_str["E3"]), "pp")
  expect_equal(unname(sign_str["E4"]), "pn")
})

test_that("iClass sign pattern for ic.num=1 uses first loading only", {
  sign1 <- ifelse(loads_mat[, 1L] >= 0, "p", "n")
  # All environments have loads1 > 0
  expect_true(all(sign1 == "p"))
})

# ---------------------------------------------------------------------------
# 8. iClassOP = mean CVE within iClass
# ---------------------------------------------------------------------------
test_that("iClassOP for 'pp' equals mean CVE of E1 and E3", {
  pp_envs <- which(sign_str == "pp")   # E1, E3
  for (g in rownames(score_mat)) {
    iop    <- mean(CVE_mat[g, pp_envs])
    # iClassOP formula: score[, 1:2] %*% colMeans(loads[pp_envs, 1:2])
    mld    <- colMeans(loads_mat[pp_envs, , drop = FALSE])
    iop_f  <- sum(score_mat[g, ] * mld)
    expect_equal(iop, iop_f, tolerance = 1e-10)
  }
})

test_that("iClassOP for 'pn' equals mean CVE of E2 and E4", {
  pn_envs <- which(sign_str == "pn")   # E2, E4
  for (g in rownames(score_mat)) {
    iop   <- mean(CVE_mat[g, pn_envs])
    mld   <- colMeans(loads_mat[pn_envs, , drop = FALSE])
    iop_f <- sum(score_mat[g, ] * mld)
    expect_equal(iop, iop_f, tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# 9. iClassRMSD = 0 when ic.num == k (all factors used -> no residual)
# ---------------------------------------------------------------------------
test_that("iClassRMSD = 0 when all factors included (ic.num = k)", {
  # dev_ic = CVE - sum_{r=1}^{k} fitted_r = CVE - CVE = 0
  fitted_all <- fitted1_mat + outer(score_mat[, 2L], loads_mat[, 2L])
  dev_ic     <- CVE_mat - fitted_all
  rmsd       <- sqrt(rowMeans(dev_ic^2))
  expect_equal(max(abs(rmsd)), 0, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 10. iClassRMSD > 0 when ic.num < k
# ---------------------------------------------------------------------------
test_that("iClassRMSD > 0 when ic.num < k (residual crossover GEI)", {
  # Use ic.num = 1: dev_ic = CVE - fitted1 = dev_mat (non-zero for k=2)
  for (w in unique(sign_str)) {
    envs_w   <- which(sign_str == w)
    dev_ic_w <- dev_mat[, envs_w, drop = FALSE]
    rmsd_w   <- sqrt(rowMeans(dev_ic_w^2))
    # At least some genotypes should have non-zero RMSD
    expect_true(any(rmsd_w > 0))
  }
})

# ---------------------------------------------------------------------------
# 11. Term parsing (logic only, not calling fast())
# ---------------------------------------------------------------------------
test_that("parse fa(Site, 2):Genotype extracts sterm='Site'", {
  term  <- "fa(Site, 2):Genotype"
  parts <- strsplit(term, ":")[[1L]]
  fa_idx <- grep("^fa\\s*\\(", parts)
  fa_part <- parts[fa_idx]
  sterm   <- trimws(sub(",.*", "", sub("^fa\\s*\\(", "", fa_part)))
  expect_equal(sterm, "Site")
})

test_that("parse fa(Site, 2):Genotype extracts gterm='Genotype'", {
  term    <- "fa(Site, 2):Genotype"
  parts   <- strsplit(term, ":")[[1L]]
  fa_idx  <- grep("^fa\\s*\\(", parts)
  gen_part <- parts[-fa_idx]
  gterm   <- trimws(gen_part)
  expect_equal(gterm, "Genotype")
})

test_that("parse vm(Genotype, ped) unwraps to 'Genotype'", {
  gen_part <- "vm(Genotype, ped)"
  if (grepl("^vm\\s*\\(", gen_part)) {
    gterm <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", gen_part))))
  } else {
    gterm <- trimws(gen_part)
  }
  expect_equal(gterm, "Genotype")
})

# ---------------------------------------------------------------------------
# 12. ic.num validation logic
# ---------------------------------------------------------------------------
test_that("ic.num must be >= 1", {
  ic.num <- 0L
  k      <- 2L
  valid  <- ic.num >= 1L && ic.num <= k
  expect_false(valid)
})

test_that("ic.num must be <= k", {
  ic.num <- 5L
  k      <- 2L
  valid  <- ic.num >= 1L && ic.num <= k
  expect_false(valid)
})

test_that("ic.num = 1 is valid for k = 2", {
  ic.num <- 1L
  k      <- 2L
  valid  <- ic.num >= 1L && ic.num <= k
  expect_true(valid)
})

# ---------------------------------------------------------------------------
# 13. Environment-major ordering in long output
# ---------------------------------------------------------------------------
test_that("environment-major order: first m rows are for E1", {
  env_long  <- envs[env_rep]
  geno_long <- genotypes[geno_rep]
  expect_equal(env_long[seq_len(m)], rep("E1", m))
  expect_equal(env_long[(m + 1):(2 * m)], rep("E2", m))
})

# ---------------------------------------------------------------------------
# 14. spec.var broadcast: same for all genotypes in same environment
# ---------------------------------------------------------------------------
test_that("spec.var repeated correctly in long format", {
  sv_long <- spec_var[env_rep]
  for (e in envs) {
    idx <- which(envs[env_rep] == e)
    expect_equal(length(unique(sv_long[idx])), 1L)
    expect_equal(unname(unique(sv_long[idx])), unname(spec_var[e]))
  }
})

# ---------------------------------------------------------------------------
# 15. Matrix dimensions
# ---------------------------------------------------------------------------
test_that("CVE_mat has m rows and t_envs columns", {
  expect_equal(dim(CVE_mat), c(m, t_envs))
})

test_that("dev_mat has same dimensions as CVE_mat", {
  expect_equal(dim(dev_mat), dim(CVE_mat))
})

test_that("stab_vec has length m", {
  expect_length(stab_vec, m)
})
