# tests/testthat/test-padTrial.R
# Tests for padTrial() — pure R, no ASReml needed.

# ---------------------------------------------------------------------------
# Shared test helpers
# ---------------------------------------------------------------------------

# Build a minimal rectangular trial layout.
#   Rows 1-4, Cols 1-5 — all cells present, single block.
make_trial <- function(nrow = 4L, ncol = 5L,
                       block = "B1", type = "DH",
                       row_col = "Row:Column") {
  pat  <- unlist(strsplit(row_col, ":", fixed = TRUE))
  grid <- expand.grid(Row    = seq_len(nrow),
                      Column = seq_len(ncol),
                      KEEP.OUT.ATTRS = FALSE)
  names(grid) <- pat
  grid$Type   <- type
  grid$Block  <- block
  grid$Geno   <- paste0("G", sprintf("%02d", seq_len(nrow(grid))))
  grid
}

# Build a layout that has guard rows (Type "Guard") around a DH core.
#   DH core: rows 2-4, cols 2-4 (3×3).
#   Extra guard rows/cols: row 1, row 5, col 1, col 5.
make_guarded_trial <- function() {
  rows <- 1L:5L
  cols <- 1L:5L
  grid <- expand.grid(Row = rows, Column = cols,
                      KEEP.OUT.ATTRS = FALSE)
  grid$Type  <- "Guard"
  # DH core
  core <- grid$Row >= 2L & grid$Row <= 4L &
          grid$Column >= 2L & grid$Column <= 4L
  grid$Type[core] <- "DH"
  grid$Block  <- "B1"
  grid$Geno   <- paste0("G", sprintf("%02d", seq_len(nrow(grid))))
  grid
}

# Build a multi-block layout (two blocks, same grid shape).
make_multi_block <- function(blocks = c("B1","B2")) {
  do.call(rbind, lapply(blocks, function(b) {
    d <- make_trial(nrow = 3L, ncol = 4L, block = b)
    d
  }))
}

# Build a layout with a missing interior cell (for padding tests).
make_missing_cell_trial <- function() {
  d <- make_trial(nrow = 3L, ncol = 3L)
  # Remove the centre cell (Row=2, Col=2) to create a gap
  d[!(d$Row == 2L & d$Column == 2L), ]
}

# ---------------------------------------------------------------------------
# 1. Input validation — pattern
# ---------------------------------------------------------------------------
test_that("invalid pattern (no colon) stops with informative error", {
  d <- make_trial()
  expect_error(padTrial(d, pattern = "RowColumn"),
               "'pattern' must be two column names separated by ':'")
})

test_that("pattern with three parts stops", {
  d <- make_trial()
  expect_error(padTrial(d, pattern = "Row:Column:Extra"),
               "'pattern' must be two column names separated by ':'")
})

test_that("pattern column not present in data stops", {
  d <- make_trial()
  expect_error(padTrial(d, pattern = "Row:XCoord"),
               "Column\\(s\\) from 'pattern' not found in data")
})

test_that("both pattern columns absent stops", {
  d <- make_trial()
  expect_error(padTrial(d, pattern = "X:Y"),
               "Column\\(s\\) from 'pattern' not found in data")
})

# ---------------------------------------------------------------------------
# 2. Input validation — split / keep
# ---------------------------------------------------------------------------
test_that("split column not in data stops", {
  d <- make_trial()
  expect_error(padTrial(d, split = "Site"),
               "'split' column\\(s\\) not found in data")
})

test_that("keep column not in data stops", {
  d <- make_trial()
  expect_error(padTrial(d, split = "Block", keep = "BadCol"),
               "'keep' column\\(s\\) not found in data")
})

# ---------------------------------------------------------------------------
# 3. Input validation — fill_value
# ---------------------------------------------------------------------------
test_that("fill_value must be a single character string", {
  d <- make_trial()
  expect_error(padTrial(d, fill_value = 123L),
               "'fill_value' must be a single character string")
})

test_that("fill_value of length > 1 stops", {
  d <- make_trial()
  expect_error(padTrial(d, fill_value = c("A","B")),
               "'fill_value' must be a single character string")
})

# ---------------------------------------------------------------------------
# 4. Input validation — type_col
# ---------------------------------------------------------------------------
test_that("type_col not in data stops", {
  d <- make_trial()
  expect_error(padTrial(d, type_col = "PlotKind"),
               "type_col.*column.*not found in data")
})

test_that("match values absent from type_col stops", {
  d <- make_trial()   # Type column contains "DH"
  expect_error(padTrial(d, match = "Line"),
               "None of the 'match' value\\(s\\) found in column")
})

# ---------------------------------------------------------------------------
# 5. Output structure
# ---------------------------------------------------------------------------
test_that("result is a data frame", {
  d <- make_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_s3_class(res, "data.frame")
})

test_that("result has an 'add' column", {
  d <- make_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_true("add" %in% names(res))
})

test_that("result has same non-add columns as input", {
  d <- make_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expected_cols <- c(names(d), "add")
  expect_true(all(expected_cols %in% names(res)))
})

test_that("'add' column only contains 'old' and 'new'", {
  d <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_true(all(res$add %in% c("old","new")))
})

# ---------------------------------------------------------------------------
# 6. No padding needed — complete grid returns all 'old'
# ---------------------------------------------------------------------------
test_that("complete grid produces no 'new' rows", {
  d <- make_trial(nrow = 3L, ncol = 4L)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_true(all(res$add == "old"))
})

test_that("complete grid: nrow(result) == nrow(input)", {
  d <- make_trial(nrow = 3L, ncol = 4L)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_equal(nrow(res), nrow(d))
})

# ---------------------------------------------------------------------------
# 7. Bounding-box subsetting (guards excluded)
# ---------------------------------------------------------------------------
test_that("guard rows outside DH bounding box are dropped", {
  d   <- make_guarded_trial()
  res <- padTrial(d, match = "DH", split = "Block", pad = FALSE,
                  verbose = FALSE)
  expect_true(all(res$Row %in% 2L:4L))
  expect_true(all(res$Column %in% 2L:4L))
})

test_that("guard plots with type 'Guard' are excluded from result", {
  d   <- make_guarded_trial()
  res <- padTrial(d, match = "DH", split = "Block", pad = FALSE,
                  verbose = FALSE)
  # Only DH plots inside bounding box remain (3×3 = 9 rows)
  expect_equal(nrow(res), 9L)
})

test_that("non-DH check plots inside bounding box are retained", {
  d      <- make_guarded_trial()
  # Place a Check in an interior cell (row 3, col 3) — inside DH bbox
  idx    <- which(d$Row == 3L & d$Column == 3L)
  d$Type[idx] <- "Check"
  res <- padTrial(d, match = "DH", split = "Block", pad = TRUE,
                  verbose = FALSE)
  expect_true("Check" %in% res$Type)
})

# ---------------------------------------------------------------------------
# 8. Padding — missing cell inserted
# ---------------------------------------------------------------------------
test_that("missing interior cell produces one 'new' row", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_equal(sum(res$add == "new"), 1L)
})

test_that("padded row has correct Row and Column coordinates", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_equal(as.integer(as.character(new_row$Row)),    2L)
  expect_equal(as.integer(as.character(new_row$Column)), 2L)
})

test_that("padded row Geno column gets fill_value", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", fill_value = "Blank", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_equal(new_row$Geno, "Blank")
})

test_that("padded row Type column gets fill_value", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", fill_value = "Blank", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_equal(new_row$Type, "Blank")
})

test_that("padded row numeric columns are NA", {
  d        <- make_missing_cell_trial()
  d$Yield  <- rnorm(nrow(d))
  res      <- padTrial(d, split = "Block", verbose = FALSE)
  new_row  <- res[res$add == "new", ]
  expect_true(is.na(new_row$Yield))
})

# ---------------------------------------------------------------------------
# 9. pad = FALSE — no rows inserted
# ---------------------------------------------------------------------------
test_that("pad = FALSE never inserts 'new' rows", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", pad = FALSE, verbose = FALSE)
  expect_true(all(res$add == "old"))
})

test_that("pad = FALSE result is a strict subset of input rows", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", pad = FALSE, verbose = FALSE)
  expect_lte(nrow(res), nrow(d))
})

# ---------------------------------------------------------------------------
# 10. split = NULL — entire dataset as one group
# ---------------------------------------------------------------------------
test_that("split = NULL works and returns a data frame", {
  d   <- make_trial()
  res <- padTrial(d, split = NULL, verbose = FALSE)
  expect_s3_class(res, "data.frame")
})

test_that("split = NULL: no extra grouping column in result", {
  d   <- make_trial()
  res <- padTrial(d, split = NULL, verbose = FALSE)
  expect_false(".padTrial_grp" %in% names(res))
})

test_that("split = NULL on complete grid: all rows marked 'old'", {
  d   <- make_trial()
  res <- padTrial(d, split = NULL, verbose = FALSE)
  expect_true(all(res$add == "old"))
})

test_that("split = NULL on missing-cell trial pads correctly", {
  d   <- make_missing_cell_trial()
  d$Block <- NULL           # remove block so we can test with no split
  res <- padTrial(d, split = NULL, verbose = FALSE)
  expect_equal(sum(res$add == "new"), 1L)
})

# ---------------------------------------------------------------------------
# 11. Single-column split
# ---------------------------------------------------------------------------
test_that("single split processes each block independently", {
  d   <- make_multi_block(blocks = c("B1","B2"))
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_true(all(c("B1","B2") %in% res$Block))
})

test_that("single split: temporary grouping key not present in result", {
  d   <- make_multi_block()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_false(".padTrial_grp" %in% names(res))
})

test_that("single split with gaps in one block pads only that block", {
  d1 <- make_trial(nrow = 3L, ncol = 3L, block = "B1")  # complete
  d2 <- make_missing_cell_trial()                         # has gap
  d2$Block <- "B2"
  d  <- rbind(d1, d2)

  res <- padTrial(d, split = "Block", verbose = FALSE)
  b1_new <- sum(res$add == "new" & res$Block == "B1")
  b2_new <- sum(res$add == "new" & res$Block == "B2")
  expect_equal(b1_new, 0L)
  expect_equal(b2_new, 1L)
})

# ---------------------------------------------------------------------------
# 12. Multi-column split
# ---------------------------------------------------------------------------
test_that("multi-column split produces correct groups", {
  d <- rbind(
    cbind(make_trial(nrow = 2L, ncol = 3L, block = "B1"), Site = "S1"),
    cbind(make_trial(nrow = 2L, ncol = 3L, block = "B1"), Site = "S2")
  )
  res <- padTrial(d, split = c("Site","Block"), keep = c("Site","Block"),
                  verbose = FALSE)
  expect_true(all(c("S1","S2") %in% res$Site))
})

test_that("multi-column split: composite key column dropped from result", {
  d <- rbind(
    cbind(make_trial(nrow = 2L, ncol = 3L, block = "B1"), Site = "S1"),
    cbind(make_trial(nrow = 2L, ncol = 3L, block = "B2"), Site = "S2")
  )
  res <- padTrial(d, split = c("Site","Block"), keep = c("Site","Block"),
                  verbose = FALSE)
  expect_false(".padTrial_grp" %in% names(res))
})

test_that("multi-column split pads each sub-group independently", {
  d1 <- make_trial(nrow = 3L, ncol = 3L, block = "B1")
  d1$Site <- "S1"
  d2 <- make_missing_cell_trial()
  d2$Block <- "B1"
  d2$Site  <- "S2"
  d  <- rbind(d1, d2)

  res <- padTrial(d, split = c("Site","Block"), keep = c("Site","Block"),
                  verbose = FALSE)
  new_s1 <- sum(res$add == "new" & res$Site == "S1")
  new_s2 <- sum(res$add == "new" & res$Site == "S2")
  expect_equal(new_s1, 0L)
  expect_equal(new_s2, 1L)
})

# ---------------------------------------------------------------------------
# 13. keep columns carried into padded rows
# ---------------------------------------------------------------------------
test_that("keep column value is preserved in padded rows", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", keep = "Block", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_equal(as.character(new_row$Block), "B1")
})

test_that("padded row 'Block' is not fill_value", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", keep = "Block",
                  fill_value = "Blank", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_true(as.character(new_row$Block) != "Blank")
})

# ---------------------------------------------------------------------------
# 14. Custom fill_value
# ---------------------------------------------------------------------------
test_that("custom fill_value appears in padded character columns", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", fill_value = "EMPTY", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  char_cols <- setdiff(names(new_row)[vapply(new_row, is.character, logical(1L))],
                       c("Row","Column","Block","add"))
  if (length(char_cols) > 0L)
    expect_true(all(new_row[, char_cols] == "EMPTY"))
})

test_that("fill_value applied to factor column adds new level", {
  d      <- make_missing_cell_trial()
  d$Type <- factor(d$Type)   # make Type a factor
  res    <- padTrial(d, split = "Block", fill_value = "EMPTY", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_true("EMPTY" %in% levels(res$Type))
  expect_equal(as.character(new_row$Type), "EMPTY")
})

# ---------------------------------------------------------------------------
# 15. Custom type_col
# ---------------------------------------------------------------------------
test_that("custom type_col name works correctly", {
  d           <- make_trial()
  names(d)[names(d) == "Type"] <- "PlotKind"
  res <- padTrial(d, type_col = "PlotKind", match = "DH",
                  split = "Block", verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true("PlotKind" %in% names(res))
})

# ---------------------------------------------------------------------------
# 16. Duplicate Row × Column detection
# ---------------------------------------------------------------------------
test_that("duplicate Row x Column in same group stops with informative error", {
  d   <- make_trial(nrow = 3L, ncol = 3L)
  dup <- d[1L, ]               # copy of first row (same Row & Column)
  d   <- rbind(d, dup)
  expect_error(padTrial(d, split = "Block", verbose = FALSE),
               "duplicate Row.*Column position")
})

test_that("duplicates across different groups do not trigger error", {
  # Same Row/Column values but different blocks — should be fine
  d1 <- make_trial(nrow = 2L, ncol = 2L, block = "B1")
  d2 <- make_trial(nrow = 2L, ncol = 2L, block = "B2")
  d  <- rbind(d1, d2)
  expect_no_error(padTrial(d, split = "Block", verbose = FALSE))
})

# ---------------------------------------------------------------------------
# 17. Output ordering
# ---------------------------------------------------------------------------
test_that("result is ordered by Row then Column ascending", {
  d   <- make_trial(nrow = 4L, ncol = 5L)
  # Shuffle input deliberately
  d   <- d[sample(nrow(d)), ]
  res <- padTrial(d, split = "Block", verbose = FALSE)
  row_num <- as.numeric(as.character(res$Row))
  col_num <- as.numeric(as.character(res$Column))
  ord     <- order(row_num, col_num)
  expect_equal(seq_len(nrow(res)), ord)
})

test_that("result ordering preserved across blocks (multi-split)", {
  d <- make_multi_block()
  d <- d[sample(nrow(d)), ]   # shuffle
  res <- padTrial(d, split = "Block", verbose = FALSE)
  row_num <- as.numeric(as.character(res$Row))
  col_num <- as.numeric(as.character(res$Column))
  ord     <- order(row_num, col_num)
  expect_equal(seq_len(nrow(res)), ord)
})

# ---------------------------------------------------------------------------
# 18. Factor re-levelling of spatial coordinates
# ---------------------------------------------------------------------------
test_that("Row in result is a factor", {
  d   <- make_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_s3_class(res$Row, "factor")
})

test_that("Column in result is a factor", {
  d   <- make_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_s3_class(res$Column, "factor")
})

test_that("Row factor levels are ascending numeric", {
  d   <- make_trial(nrow = 4L, ncol = 3L)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  lvls <- as.integer(levels(res$Row))
  expect_equal(lvls, sort(lvls))
})

test_that("Column factor levels are ascending numeric", {
  d   <- make_trial(nrow = 3L, ncol = 5L)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  lvls <- as.integer(levels(res$Column))
  expect_equal(lvls, sort(lvls))
})

test_that("padded row coordinates appear as valid factor levels", {
  d   <- make_missing_cell_trial()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_true(as.character(new_row$Row)    %in% levels(res$Row))
  expect_true(as.character(new_row$Column) %in% levels(res$Column))
})

# ---------------------------------------------------------------------------
# 19. Multiple missing cells
# ---------------------------------------------------------------------------
test_that("multiple missing cells all padded correctly", {
  # 4×4 grid — remove 3 cells
  d <- make_trial(nrow = 4L, ncol = 4L)
  remove_idx <- which((d$Row == 1L & d$Column == 3L) |
                        (d$Row == 2L & d$Column == 2L) |
                        (d$Row == 4L & d$Column == 4L))
  d <- d[-remove_idx, ]
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_equal(sum(res$add == "new"), 3L)
  expect_equal(nrow(res), 4L * 4L)
})

# ---------------------------------------------------------------------------
# 20. Multiple match types
# ---------------------------------------------------------------------------
test_that("multiple match values use the union bounding box", {
  d <- make_trial(nrow = 5L, ncol = 5L)
  d$Type[d$Row == 1L] <- "Check"
  # DH rows: 2-5; Check row 1 is also 'match'-able
  res <- padTrial(d, match = c("DH","Check"),
                  split = "Block", pad = FALSE, verbose = FALSE)
  # Bounding box covers rows 1-5 since Check row 1 is in match
  expect_true(1L %in% as.integer(as.character(res$Row)))
})

# ---------------------------------------------------------------------------
# 21. verbose output (message emitted)
# ---------------------------------------------------------------------------
test_that("verbose = TRUE emits a message", {
  d <- make_trial()
  expect_message(padTrial(d, split = "Block", verbose = TRUE),
                 regexp = NULL)
})

test_that("verbose = FALSE emits no messages", {
  d <- make_trial()
  expect_no_message(padTrial(d, split = "Block", verbose = FALSE))
})

# ---------------------------------------------------------------------------
# 22. Group with no matching rows returned unchanged (verbose warning path)
# ---------------------------------------------------------------------------
test_that("group with no 'match' rows returned unchanged without error", {
  d1 <- make_trial(nrow = 2L, ncol = 2L, block = "B1", type = "DH")
  d2 <- make_trial(nrow = 2L, ncol = 2L, block = "B2", type = "Check")
  d  <- rbind(d1, d2)

  # Only DH is in 'match'; group B2 has no DH rows — should not error
  expect_no_error(
    suppressMessages(
      padTrial(d, match = "DH", split = "Block", verbose = TRUE)
    )
  )
})

test_that("group with no 'match' rows still has 'add' column in result", {
  d1 <- make_trial(nrow = 2L, ncol = 2L, block = "B1", type = "DH")
  d2 <- make_trial(nrow = 2L, ncol = 2L, block = "B2", type = "Check")
  d  <- rbind(d1, d2)
  res <- suppressMessages(
    padTrial(d, match = "DH", split = "Block", verbose = TRUE)
  )
  expect_true("add" %in% names(res))
})

# ---------------------------------------------------------------------------
# 23. keep = split default behaviour
# ---------------------------------------------------------------------------
test_that("default keep equals split: Block in padded rows when split = 'Block'", {
  d   <- make_missing_cell_trial()
  # Default keep = split = "Block"
  res <- padTrial(d, split = "Block", verbose = FALSE)
  new_row <- res[res$add == "new", ]
  expect_equal(as.character(new_row$Block), "B1")
})

# ---------------------------------------------------------------------------
# 24. Row count sanity checks
# ---------------------------------------------------------------------------
test_that("result nrow equals n_old + n_padded", {
  d       <- make_missing_cell_trial()
  n_old   <- nrow(d)
  res     <- padTrial(d, split = "Block", verbose = FALSE)
  n_new   <- sum(res$add == "new")
  expect_equal(nrow(res), n_old + n_new)
})

test_that("complete 3×3 grid: result has exactly 9 rows", {
  d   <- make_trial(nrow = 3L, ncol = 3L)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_equal(nrow(res), 9L)
})

# ---------------------------------------------------------------------------
# 25. Non-sequential (non-integer) coordinate values
# ---------------------------------------------------------------------------
test_that("non-sequential row numbers handled correctly", {
  d        <- make_trial(nrow = 4L, ncol = 3L)
  d$Row    <- d$Row * 2L    # rows: 2, 4, 6, 8  (gaps but consistent)
  d$Column <- d$Column * 3L # cols: 3, 6, 9
  # No interior gaps within the bounding box, so no padding needed
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_true(all(res$add == "old"))
})

test_that("non-sequential coords: padding is table-based on observed levels", {
  # Rows 1 and 3 only — row 2 never appears so the table has no 'row 2' level.
  # padTrial uses table() which only sees observed coordinate values, so no
  # gap-filling between 1 and 3.  All 4 rows should be marked 'old'.
  d <- data.frame(Row = c(1L,1L,3L,3L),
                  Column = c(1L,2L,1L,2L),
                  Type = "DH", Block = "B1", Geno = paste0("G",1:4),
                  stringsAsFactors = FALSE)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_true(all(res$add == "old"))
})

# ---------------------------------------------------------------------------
# 26. Factor input columns preserved
# ---------------------------------------------------------------------------
test_that("factor Type column preserved in old rows", {
  d      <- make_trial()
  d$Type <- factor(d$Type)
  res    <- padTrial(d, split = "Block", verbose = FALSE)
  expect_s3_class(res$Type, "factor")
})

test_that("factor Block column preserved in old rows", {
  d       <- make_trial()
  d$Block <- factor(d$Block)
  res     <- padTrial(d, split = "Block", verbose = FALSE)
  expect_s3_class(res$Block, "factor")
})

# ---------------------------------------------------------------------------
# 27. Bounding box when all plots are match type (no guard trimming)
# ---------------------------------------------------------------------------
test_that("all-DH data: bounding box covers entire grid, nrow unchanged", {
  d   <- make_trial(nrow = 4L, ncol = 4L, type = "DH")
  res <- padTrial(d, match = "DH", split = "Block", verbose = FALSE)
  expect_equal(nrow(res), 4L * 4L)
})

# ---------------------------------------------------------------------------
# 28. Large bounding box with multiple padded cells (multi-block)
# ---------------------------------------------------------------------------
test_that("two-block trial with gaps padded to full rectangles", {
  mk_with_gap <- function(blk) {
    d <- make_trial(nrow = 4L, ncol = 4L, block = blk)
    # Remove 2 cells per block
    d[!(d$Row == 2L & d$Column == 3L) &
      !(d$Row == 3L & d$Column == 1L), ]
  }
  d   <- rbind(mk_with_gap("B1"), mk_with_gap("B2"))
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_equal(sum(res$add == "new"), 4L)   # 2 gaps × 2 blocks
  expect_equal(nrow(res), 4L * 4L * 2L)
})
