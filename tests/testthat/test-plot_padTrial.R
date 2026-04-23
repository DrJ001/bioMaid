# tests/testthat/test-plot_padTrial.R
# Tests for plot_padTrial() — pure R, no ASReml needed.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Minimal complete 3×3 trial — no missing cells, one block
make_trial_data <- function(nrow = 3L, ncol = 3L, block = "B1",
                             type = "DH", seed = 1L) {
  set.seed(seed)
  d <- expand.grid(Row    = seq_len(nrow),
                   Column = seq_len(ncol),
                   KEEP.OUT.ATTRS = FALSE)
  d$Type  <- type
  d$Block <- block
  d$Geno  <- paste0("G", sprintf("%02d", seq_len(nrow(d))))
  d
}

# Trial with one missing interior cell (Row=2, Col=2)
make_missing_data <- function(block = "B1") {
  d <- make_trial_data(nrow = 3L, ncol = 3L, block = block)
  d[!(d$Row == 2L & d$Column == 2L), ]
}

# Trial with guard rows surrounding a DH core
make_guarded_data <- function() {
  rows <- 1L:5L; cols <- 1L:5L
  d <- expand.grid(Row = rows, Column = cols, KEEP.OUT.ATTRS = FALSE)
  d$Type  <- "Guard"
  core    <- d$Row >= 2L & d$Row <= 4L & d$Column >= 2L & d$Column <= 4L
  d$Type[core] <- "DH"
  d$Block <- "B1"
  d$Geno  <- paste0("G", sprintf("%02d", seq_len(nrow(d))))
  d
}

# Run padTrial and return result
make_pt_result <- function(data = NULL, split = "Block",
                            pad = TRUE, verbose = FALSE) {
  if (is.null(data)) data <- make_missing_data()
  padTrial(data, split = split, pad = pad, verbose = verbose)
}

# ---------------------------------------------------------------------------
# 1. Basic return types
# ---------------------------------------------------------------------------
test_that("plot_padTrial returns a ggplot object", {
  res <- make_pt_result()
  expect_s3_class(plot_padTrial(res), "ggplot")
})

test_that("return_data = TRUE returns a data frame", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_s3_class(df, "data.frame")
})

test_that("plot_padTrial is returned invisibly", {
  res <- make_pt_result()
  expect_invisible(plot_padTrial(res))
})

# ---------------------------------------------------------------------------
# 2. Data frame structure (return_data = TRUE)
# ---------------------------------------------------------------------------
test_that("return_data has required columns", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expected <- c("row_num", "col_num", "fill_group", "panel",
                "group", "is_new", "label_text")
  expect_true(all(expected %in% names(df)))
})

test_that("panel column has levels Before and After", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_equal(levels(df$panel), c("Before", "After"))
})

test_that("both Before and After panels are present", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_true("Before" %in% df$panel)
  expect_true("After"  %in% df$panel)
})

test_that("row_num and col_num are numeric", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_type(df$row_num, "double")
  expect_type(df$col_num, "double")
})

test_that("is_new is logical", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_type(df$is_new, "logical")
})

# ---------------------------------------------------------------------------
# 3. Padded rows appear correctly in After panel
# ---------------------------------------------------------------------------
test_that("padded rows (add=new) appear as 'Padded' in After panel", {
  res <- make_pt_result()   # has one missing cell -> one new row
  df  <- plot_padTrial(res, return_data = TRUE)
  aft <- df[df$panel == "After", ]
  expect_true("Padded" %in% aft$fill_group)
})

test_that("padded rows have is_new = TRUE", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  padded <- df[df$fill_group == "Padded", ]
  expect_true(all(padded$is_new))
})

test_that("padded rows only appear in the After panel", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  bef_padded <- df[df$panel == "Before" & df$fill_group == "Padded", ]
  expect_equal(nrow(bef_padded), 0L)
})

test_that("number of Padded rows equals number of new rows in result", {
  res <- make_pt_result()
  n_new <- sum(res$add == "new")
  df    <- plot_padTrial(res, return_data = TRUE)
  n_pad <- sum(df$fill_group == "Padded")
  expect_equal(n_pad, n_new)
})

# ---------------------------------------------------------------------------
# 4. data = NULL (Before reconstructed from old rows)
# ---------------------------------------------------------------------------
test_that("data=NULL: Before panel contains only old rows", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  bef <- df[df$panel == "Before", ]
  expect_true(all(!bef$is_new))
  expect_false("Padded" %in% bef$fill_group)
})

test_that("data=NULL: Before row count equals number of old rows in result", {
  res   <- make_pt_result()
  n_old <- sum(res$add == "old")
  df    <- plot_padTrial(res, return_data = TRUE)
  n_bef <- sum(df$panel == "Before")
  expect_equal(n_bef, n_old)
})

test_that("data=NULL: After row count equals nrow(result)", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_equal(sum(df$panel == "After"), nrow(res))
})

# ---------------------------------------------------------------------------
# 5. data supplied (Before from original data)
# ---------------------------------------------------------------------------
test_that("data supplied: Before panel row count equals nrow(data)", {
  orig <- make_missing_data()
  res  <- make_pt_result(data = orig)
  df   <- plot_padTrial(res, data = orig, return_data = TRUE)
  expect_equal(sum(df$panel == "Before"), nrow(orig))
})

test_that("data supplied: Before panel covers all original row/col coords", {
  orig <- make_missing_data()
  res  <- make_pt_result(data = orig)
  df   <- plot_padTrial(res, data = orig, return_data = TRUE)
  bef  <- df[df$panel == "Before", ]
  expect_equal(sort(unique(bef$row_num)),
               sort(unique(as.numeric(orig$Row))))
  expect_equal(sort(unique(bef$col_num)),
               sort(unique(as.numeric(orig$Column))))
})

test_that("data supplied with guard rows: Before includes guard row coords", {
  orig <- make_guarded_data()
  res  <- padTrial(orig, match = "DH", split = "Block",
                   verbose = FALSE)
  df   <- plot_padTrial(res, data = orig, return_data = TRUE)
  bef  <- df[df$panel == "Before", ]
  # Guard rows include Row=1 and Row=5 which are outside the DH bounding box
  expect_true(1L %in% bef$row_num)
  expect_true(5L %in% bef$row_num)
})

test_that("data supplied with guard rows: After excludes guard row coords", {
  orig <- make_guarded_data()
  res  <- padTrial(orig, match = "DH", split = "Block",
                   verbose = FALSE)
  df   <- plot_padTrial(res, data = orig, return_data = TRUE)
  aft  <- df[df$panel == "After", ]
  # After panel is bounded by DH bounding box (rows 2-4, cols 2-4)
  expect_false(1L %in% aft$row_num)
  expect_false(5L %in% aft$row_num)
})

# ---------------------------------------------------------------------------
# 6. split argument
# ---------------------------------------------------------------------------
test_that("split = 'Block' produces group column values", {
  d1  <- make_missing_data("B1")
  d2  <- make_missing_data("B2")
  d2$Geno <- paste0("H", sprintf("%02d", seq_len(nrow(d2))))
  d   <- rbind(d1, d2)
  res <- padTrial(d, split = "Block", verbose = FALSE)
  df  <- plot_padTrial(res, split = "Block", return_data = TRUE)
  expect_true(all(c("B1", "B2") %in% df$group))
})

test_that("no split: group column is 'All'", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_true(all(df$group == "All"))
})

test_that("multi-column split produces middle-dot composite group labels", {
  orig <- make_missing_data()
  orig$Site <- "S1"
  res <- padTrial(orig, split = c("Site", "Block"),
                  keep = c("Site", "Block"), verbose = FALSE)
  df  <- plot_padTrial(res, split = c("Site", "Block"),
                       return_data = TRUE)
  expect_true(any(grepl("\u00b7", df$group)))
})

# ---------------------------------------------------------------------------
# 7. label argument
# ---------------------------------------------------------------------------
test_that("label = NULL produces empty label_text", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, label = NULL, return_data = TRUE)
  expect_true(all(df$label_text == ""))
})

test_that("label = 'Geno' populates label_text for old rows", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, label = "Geno", return_data = TRUE)
  aft_old <- df[df$panel == "After" & !df$is_new, ]
  expect_true(any(nchar(aft_old$label_text) > 0L))
})

test_that("label = 'Geno': padded tiles have empty label_text", {
  res <- make_pt_result()
  df  <- plot_padTrial(res, label = "Geno", return_data = TRUE)
  padded <- df[df$is_new, ]
  expect_true(all(padded$label_text == ""))
})

# ---------------------------------------------------------------------------
# 8. Complete grid (no new rows) still plots
# ---------------------------------------------------------------------------
test_that("complete grid (no padding) plots without error", {
  d   <- make_trial_data()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  expect_s3_class(plot_padTrial(res), "ggplot")
})

test_that("complete grid: no Padded fill_group in data", {
  d   <- make_trial_data()
  res <- padTrial(d, split = "Block", verbose = FALSE)
  df  <- plot_padTrial(res, return_data = TRUE)
  expect_false("Padded" %in% df$fill_group)
})

# ---------------------------------------------------------------------------
# 9. Internal helper: .pt_parse_pattern
# ---------------------------------------------------------------------------
test_that(".pt_parse_pattern returns correct row_col and col_col", {
  rc <- biomAid:::.pt_parse_pattern("Row:Column")
  expect_equal(rc$row_col, "Row")
  expect_equal(rc$col_col, "Column")
})

test_that(".pt_parse_pattern handles whitespace", {
  rc <- biomAid:::.pt_parse_pattern("  Row : Column ")
  expect_equal(rc$row_col, "Row")
  expect_equal(rc$col_col, "Column")
})

test_that(".pt_parse_pattern errors on missing colon", {
  expect_error(biomAid:::.pt_parse_pattern("RowColumn"),
               "'pattern' must be two column names separated by ':'")
})

test_that(".pt_parse_pattern errors on three parts", {
  expect_error(biomAid:::.pt_parse_pattern("Row:Col:Extra"),
               "'pattern' must be two column names separated by ':'")
})

# ---------------------------------------------------------------------------
# 10. Input validation
# ---------------------------------------------------------------------------
test_that("non-data-frame result errors", {
  expect_error(plot_padTrial("not_a_df"), "'result' must be a data frame")
})

test_that("result without 'add' column errors", {
  res <- make_pt_result()
  res$add <- NULL
  expect_error(plot_padTrial(res), "'add' column")
})

test_that("result with invalid 'add' values errors", {
  res     <- make_pt_result()
  res$add <- rep("bad", nrow(res))
  expect_error(plot_padTrial(res), "\"old\" and \"new\"")
})

test_that("invalid pattern (no colon) errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, pattern = "RowColumn"),
               "'pattern' must be two column names")
})

test_that("pattern column not in result errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, pattern = "Row:XCoord"),
               "not found in 'result'")
})

test_that("type_col not in result errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, type_col = "PlotKind"),
               "'type_col' column.*not found")
})

test_that("split column not in result errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, split = "Site"),
               "'split' column\\(s\\) not found in 'result'")
})

test_that("split column not in data errors", {
  orig <- make_missing_data()
  res  <- make_pt_result(data = orig)
  expect_error(
    plot_padTrial(res, data = orig, split = "Site"),
    "'split' column\\(s\\) not found in 'result'"
  )
})

test_that("label not a string errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, label = 123),
               "'label' must be a single character string")
})

test_that("label column not in result errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, label = "NoSuchCol"),
               "'label' column.*not found")
})

test_that("non-data-frame data errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, data = "bad"),
               "'data' must be a data frame")
})

test_that("data missing pattern column errors", {
  orig      <- make_missing_data()
  bad_data  <- orig
  bad_data$Row <- NULL
  res       <- make_pt_result(data = orig)
  expect_error(plot_padTrial(res, data = bad_data),
               "not found in 'data'")
})

test_that("non-theme theme errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, theme = "bw"),
               "'theme' must be a ggplot2 theme")
})

test_that("non-logical return_data errors", {
  res <- make_pt_result()
  expect_error(plot_padTrial(res, return_data = "yes"),
               "'return_data' must be TRUE or FALSE")
})

# ---------------------------------------------------------------------------
# 11. Custom type_col name
# ---------------------------------------------------------------------------
test_that("custom type_col works correctly", {
  d          <- make_trial_data()
  names(d)[names(d) == "Type"] <- "PlotKind"
  res        <- padTrial(d, type_col = "PlotKind", split = "Block",
                          verbose = FALSE)
  df         <- plot_padTrial(res, type_col = "PlotKind",
                               return_data = TRUE)
  expect_true("DH" %in% df$fill_group)
})

# ---------------------------------------------------------------------------
# 12. Mixed plot types
# ---------------------------------------------------------------------------
test_that("multiple type values all appear in fill_group", {
  d        <- make_trial_data(nrow = 4L, ncol = 4L)
  d$Type[d$Row == 1L] <- "Check"
  res      <- padTrial(d, match = c("DH","Check"),
                       split = "Block", verbose = FALSE)
  df       <- plot_padTrial(res, return_data = TRUE)
  aft_types <- unique(df$fill_group[df$panel == "After"])
  expect_true("DH"    %in% aft_types)
  expect_true("Check" %in% aft_types)
})
