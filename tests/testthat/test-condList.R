# tests/testthat/test-condList.R
# Dedicated tests for .condList() — shared by randomRegressMV & fixedRegressMV

test_that(".condList: baseline structure", {
  cl <- biomAid:::.condList(c("N0","N1","N2","N3"), "baseline", NULL)
  expect_null(cl[["N0"]])
  expect_equal(cl[["N1"]], "N0")
  expect_equal(cl[["N2"]], "N0")
  expect_equal(cl[["N3"]], "N0")
})

test_that(".condList: sequential is Gram-Schmidt (LDL' Cholesky)", {
  cl <- biomAid:::.condList(c("A","B","C","D"), "sequential", NULL)
  expect_null(cl[["A"]])
  expect_equal(cl[["B"]], "A")
  expect_equal(cl[["C"]], c("A","B"))
  expect_equal(cl[["D"]], c("A","B","C"))
})

test_that(".condList: partial — each vs all others", {
  cl <- biomAid:::.condList(c("X","Y","Z"), "partial", NULL)
  expect_setequal(cl[["X"]], c("Y","Z"))
  expect_setequal(cl[["Y"]], c("X","Z"))
  expect_setequal(cl[["Z"]], c("X","Y"))
})

test_that(".condList: two-level baseline has one conditioned treatment", {
  cl <- biomAid:::.condList(c("T0","T1"), "baseline", NULL)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
})

test_that(".condList: custom with full specification", {
  cond <- list(T0 = NULL, T1 = "T0", T2 = c("T0","T1"))
  cl   <- biomAid:::.condList(c("T0","T1","T2"), "custom", cond)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_equal(cl[["T2"]], c("T0","T1"))
})

test_that(".condList: custom with partial specification (unspecified -> NULL)", {
  cond <- list(T1 = "T0")   # T0 and T2 not specified -> remain NULL
  cl   <- biomAid:::.condList(c("T0","T1","T2"), "custom", cond)
  expect_null(cl[["T0"]])
  expect_equal(cl[["T1"]], "T0")
  expect_null(cl[["T2"]])
})

test_that(".condList custom: missing cond errors with message", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", NULL),
    "'cond' must be supplied"
  )
})

test_that(".condList custom: unnamed cond list errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", list("A")),
    "named list"
  )
})

test_that(".condList custom: unrecognised name in cond errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", list(Z = NULL)),
    "not found in 'levs'"
  )
})

test_that(".condList custom: conditioning level not in levs errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", list(B = "C")),
    "not in 'levs'"
  )
})

test_that(".condList custom: self-referential conditioning errors", {
  expect_error(
    biomAid:::.condList(c("A","B"), "custom", list(B = c("A","B"))),
    "cannot be in its own conditioning set"
  )
})

test_that(".condList: number of conditioned treatments correct", {
  # baseline: k-1 conditioned
  cl_b <- biomAid:::.condList(c("A","B","C"), "baseline", NULL)
  nc_b <- sum(!vapply(cl_b, is.null, logical(1L)))
  expect_equal(nc_b, 2L)

  # sequential: k-1 conditioned
  cl_s <- biomAid:::.condList(c("A","B","C"), "sequential", NULL)
  nc_s <- sum(!vapply(cl_s, is.null, logical(1L)))
  expect_equal(nc_s, 2L)

  # partial: all k conditioned
  cl_p <- biomAid:::.condList(c("A","B","C"), "partial", NULL)
  nc_p <- sum(!vapply(cl_p, is.null, logical(1L)))
  expect_equal(nc_p, 3L)
})
