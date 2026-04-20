# tests/testthat/test-simTrialData.R
# Tests for simTrialData() — pure R, no ASReml needed.

# ---------------------------------------------------------------------------
# 1. Basic output structure
# ---------------------------------------------------------------------------
test_that("simTrialData returns a named list with 'data' and 'params'", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("N0","N1"), nrep = 2L,
                      seed = 1L, verbose = FALSE)
  expect_named(out, c("data","params"))
  expect_s3_class(out$data, "data.frame")
  expect_true(is.list(out$params))
})

# ---------------------------------------------------------------------------
# 2. Data dimensions
# ---------------------------------------------------------------------------
test_that("total plot count = nvar * nsite * ntreat * nrep", {
  nvar <- 4L; nsite <- 3L; ntreat <- 2L; nrep <- 2L
  out <- simTrialData(nvar = nvar, nsite = nsite,
                      treatments = c("T0","T1"), nrep = nrep,
                      seed = 2L, verbose = FALSE)
  expect_equal(nrow(out$data), nvar * nsite * ntreat * nrep)
})

test_that("data has required columns", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 3L, verbose = FALSE)
  expected_cols <- c("Site","Treatment","TSite","Variety","Rep","Row","Column","yield")
  expect_true(all(expected_cols %in% names(out$data)))
})

# ---------------------------------------------------------------------------
# 3. Factor columns
# ---------------------------------------------------------------------------
test_that("Site, Treatment, TSite, Variety, Rep are factors", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 4L, verbose = FALSE)
  d <- out$data
  expect_s3_class(d$Site,      "factor")
  expect_s3_class(d$Treatment, "factor")
  expect_s3_class(d$TSite,     "factor")
  expect_s3_class(d$Variety,   "factor")
  expect_s3_class(d$Rep,       "factor")
})

# ---------------------------------------------------------------------------
# 4. Factor levels match specified inputs
# ---------------------------------------------------------------------------
test_that("Treatment factor levels match 'treatments' argument", {
  trt <- c("Dry","Irr","Wet")
  out <- simTrialData(nvar = 6L, nsite = 2L, treatments = trt,
                      nrep = 2L, seed = 5L, verbose = FALSE)
  expect_equal(levels(out$data$Treatment), trt)
})

test_that("number of Site levels equals nsite", {
  out <- simTrialData(nvar = 4L, nsite = 5L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 6L, verbose = FALSE)
  expect_equal(nlevels(out$data$Site), 5L)
})

test_that("number of Variety levels equals nvar", {
  out <- simTrialData(nvar = 9L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 7L, verbose = FALSE)
  expect_equal(nlevels(out$data$Variety), 9L)
})

# ---------------------------------------------------------------------------
# 5. TSite: treatment-site combinations
# ---------------------------------------------------------------------------
test_that("TSite has nsite * ntreat levels", {
  nsite <- 3L; ntreat <- 3L
  out <- simTrialData(nvar = 4L, nsite = nsite,
                      treatments = c("T0","T1","T2"), nrep = 2L,
                      seed = 8L, verbose = FALSE)
  expect_equal(nlevels(out$data$TSite), nsite * ntreat)
})

test_that("TSite labels contain treatment and site joined by sep", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      sep = "-", seed = 9L, verbose = FALSE)
  tsite_levs <- levels(out$data$TSite)
  expect_true(all(grepl("-", tsite_levs, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 6. Reproducibility with seed
# ---------------------------------------------------------------------------
test_that("same seed produces identical data", {
  out1 <- simTrialData(nvar = 4L, nsite = 2L,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 42L, verbose = FALSE)
  out2 <- simTrialData(nvar = 4L, nsite = 2L,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 42L, verbose = FALSE)
  expect_identical(out1$data$yield, out2$data$yield)
})

test_that("different seeds produce different data", {
  out1 <- simTrialData(nvar = 4L, nsite = 2L,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 1L, verbose = FALSE)
  out2 <- simTrialData(nvar = 4L, nsite = 2L,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 2L, verbose = FALSE)
  expect_false(identical(out1$data$yield, out2$data$yield))
})

# ---------------------------------------------------------------------------
# 7. params structure
# ---------------------------------------------------------------------------
test_that("params has all expected elements", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 10L, verbose = FALSE)
  expected_params <- c("beta","sigr","treat_effects","site_means","g_base","g_arr")
  expect_true(all(expected_params %in% names(out$params)))
})

test_that("params$beta has nsite rows and ntreat-1 columns", {
  nsite <- 3L; ntreat <- 3L
  out <- simTrialData(nvar = 4L, nsite = nsite,
                      treatments = c("T0","T1","T2"), nrep = 2L,
                      seed = 11L, verbose = FALSE)
  expect_equal(dim(out$params$beta), c(nsite, ntreat - 1L))
})

test_that("params$beta row names are site labels", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0","T1"), nrep = 2L,
                      site_prefix = "Env",
                      seed = 12L, verbose = FALSE)
  rn <- rownames(out$params$beta)
  expect_true(all(grepl("^Env", rn)))
})

test_that("params$sigr has length ntreat - 1", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1","T2","T3"), nrep = 2L,
                      seed = 13L, verbose = FALSE)
  expect_length(out$params$sigr, 3L)
})

test_that("params$g_arr has dim nvar x ntreat x nsite", {
  nvar <- 6L; nsite <- 2L; ntreat <- 2L
  out  <- simTrialData(nvar = nvar, nsite = nsite,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 14L, verbose = FALSE)
  expect_equal(dim(out$params$g_arr), c(nvar, ntreat, nsite))
})

# ---------------------------------------------------------------------------
# 8. Error handling
# ---------------------------------------------------------------------------
test_that("< 2 treatments errors", {
  expect_error(
    simTrialData(treatments = "T0", verbose = FALSE),
    "'treatments' must have at least 2"
  )
})

test_that("sep character in treatment label errors", {
  expect_error(
    simTrialData(treatments = c("T-0","T1"), sep = "-", verbose = FALSE),
    "'sep'.*must not appear"
  )
})

test_that("rows_per_strip * cols_per_strip != nvar errors", {
  expect_error(
    simTrialData(nvar = 10L, rows_per_strip = 3L, cols_per_strip = 4L,
                 treatments = c("T0","T1"), verbose = FALSE),
    "rows_per_strip \\* cols_per_strip"
  )
})

test_that("treat_effects wrong length errors", {
  expect_error(
    simTrialData(treatments = c("T0","T1","T2"),
                 treat_effects = c(0, 100),
                 verbose = FALSE),
    "'treat_effects' must have length"
  )
})

test_that("sigr wrong length errors", {
  expect_error(
    simTrialData(treatments = c("T0","T1","T2"),
                 sigr = c(100),
                 verbose = FALSE),
    "'sigr' must have length"
  )
})

test_that("beta_min >= beta_max errors", {
  expect_error(
    simTrialData(treatments = c("T0","T1"),
                 beta_min = 1.5, beta_max = 0.5,
                 verbose = FALSE),
    "'beta_min' must be strictly less"
  )
})

# ---------------------------------------------------------------------------
# 9. Data ordering: sorted by Site, Row, Column
# ---------------------------------------------------------------------------
test_that("data is sorted by Site then Row then Column", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 20L, verbose = FALSE)
  d   <- out$data
  ord <- order(d$Site, d$Row, d$Column)
  expect_equal(seq_len(nrow(d)), ord)
})

# ---------------------------------------------------------------------------
# 10. yield is numeric and finite
# ---------------------------------------------------------------------------
test_that("yield is numeric with no NAs or infinites", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      seed = 21L, verbose = FALSE)
  expect_true(is.numeric(out$data$yield))
  expect_false(anyNA(out$data$yield))
  expect_true(all(is.finite(out$data$yield)))
})

# ---------------------------------------------------------------------------
# 11. Rep factor has nrep levels
# ---------------------------------------------------------------------------
test_that("Rep factor has nrep levels", {
  nrep <- 4L
  out  <- simTrialData(nvar = 4L, nsite = 2L,
                       treatments = c("T0","T1"), nrep = nrep,
                       seed = 22L, verbose = FALSE)
  expect_equal(nlevels(out$data$Rep), nrep)
})

# ---------------------------------------------------------------------------
# 12. Variety prefix
# ---------------------------------------------------------------------------
test_that("Variety labels start with variety_prefix", {
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      variety_prefix = "Genotype",
                      seed = 23L, verbose = FALSE)
  expect_true(all(grepl("^Genotype", levels(out$data$Variety))))
})

# ---------------------------------------------------------------------------
# 13. site_prefix
# ---------------------------------------------------------------------------
test_that("Site labels start with site_prefix", {
  out <- simTrialData(nvar = 4L, nsite = 3L,
                      treatments = c("T0","T1"), nrep = 2L,
                      site_prefix = "Trial",
                      seed = 24L, verbose = FALSE)
  expect_true(all(grepl("^Trial", levels(out$data$Site))))
})

# ---------------------------------------------------------------------------
# 14. treat_effects shift means
# ---------------------------------------------------------------------------
test_that("positive treat_effects increases mean yield for that treatment", {
  out <- simTrialData(nvar = 10L, nsite = 5L,
                      treatments = c("T0","T1"),
                      treat_effects = c(0, 1000),
                      sigma_base = 10, error_sd = 50, rep_sd = 20,
                      row_sd = 10, col_sd = 10,
                      seed = 25L, verbose = FALSE)
  means <- tapply(out$data$yield, out$data$Treatment, mean)
  expect_gt(means["T1"], means["T0"])
})

# ---------------------------------------------------------------------------
# 15. outfile writes a CSV
# ---------------------------------------------------------------------------
test_that("outfile writes CSV and file exists", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  simTrialData(nvar = 4L, nsite = 2L,
               treatments = c("T0","T1"), nrep = 2L,
               seed = 26L, outfile = tmp, verbose = FALSE)
  expect_true(file.exists(tmp))
  d <- read.csv(tmp, stringsAsFactors = FALSE)
  expect_s3_class(d, "data.frame")
  expect_true("yield" %in% names(d))
})

# ---------------------------------------------------------------------------
# 16. Prime nvar strips layout: 1 x nvar dimensions
# ---------------------------------------------------------------------------
test_that("prime nvar results in 1 x nvar strip dimensions", {
  # .best_dims internal function returns c(rows=1, cols=nvar) for prime nvar
  # We verify indirectly: Row values should max at nrep * 1 = nrep
  suppressWarnings({
    out <- simTrialData(nvar = 7L, nsite = 2L,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 27L, verbose = FALSE)
  })
  # rows_per_strip = 1 when nvar is prime, so max Row per site = nrep * 1
  site_rows <- tapply(out$data$Row, out$data$Site, max)
  expect_true(all(site_rows == 2L))   # nrep = 2, rows_per_strip = 1
})

# ---------------------------------------------------------------------------
# 17. g_base has length nvar
# ---------------------------------------------------------------------------
test_that("g_base has length nvar", {
  nvar <- 8L
  out  <- simTrialData(nvar = nvar, nsite = 2L,
                       treatments = c("T0","T1"), nrep = 2L,
                       seed = 28L, verbose = FALSE)
  expect_length(out$params$g_base, nvar)
})

# ---------------------------------------------------------------------------
# 18. Row and Column range
# ---------------------------------------------------------------------------
test_that("Row values in valid range", {
  nvar <- 4L; nsite <- 2L; nrep <- 2L
  out <- simTrialData(nvar = nvar, nsite = nsite,
                      treatments = c("T0","T1"), nrep = nrep,
                      seed = 29L, verbose = FALSE)
  d <- out$data
  # rows_per_strip * nrep = total rows per site
  dims <- c(rows = as.integer(floor(sqrt(nvar))),
            cols = as.integer(nvar %/% floor(sqrt(nvar))))
  nrow_e <- nrep * dims["rows"]
  expect_true(all(d$Row >= 1L & d$Row <= nrow_e))
})

# ---------------------------------------------------------------------------
# 19. Custom rows_per_strip and cols_per_strip
# ---------------------------------------------------------------------------
test_that("explicit rows_per_strip x cols_per_strip works when product = nvar", {
  out <- simTrialData(nvar = 6L, nsite = 2L,
                      treatments = c("T0","T1"), nrep = 2L,
                      rows_per_strip = 2L, cols_per_strip = 3L,
                      seed = 30L, verbose = FALSE)
  expect_equal(nrow(out$data), 6L * 2L * 2L * 2L)
})

# ---------------------------------------------------------------------------
# 20. params$treat_effects names match treatment vector
# ---------------------------------------------------------------------------
test_that("treat_effects names match treatment argument", {
  trt <- c("N0","N1","N2")
  out <- simTrialData(nvar = 4L, nsite = 2L,
                      treatments = trt, nrep = 2L,
                      seed = 31L, verbose = FALSE)
  expect_equal(names(out$params$treat_effects), trt)
})
