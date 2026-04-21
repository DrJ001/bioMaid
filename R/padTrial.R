#' Extract and Pad a Sub-Trial from a Field Trial Layout
#'
#' @description
#' Given a field trial data frame that mixes multiple plot types (e.g. test
#' lines, checks, guard rows), `padTrial()` extracts the rectangular sub-trial
#' occupied by a target plot type and, optionally, inserts placeholder rows for
#' any missing grid positions within that rectangle.
#'
#' **Steps performed per block:**
#' \enumerate{
#'   \item Identify the bounding box (min/max of the two spatial coordinates)
#'         for all plots matching `match`.
#'   \item Subset the data to that bounding rectangle (dropping checks and
#'         guards that fall outside it).
#'   \item If `pad = TRUE`, detect any Row × Column cells that are absent
#'         within the bounding box and insert placeholder rows for them,
#'         carrying identifier columns specified in `keep` and writing
#'         `fill_value` into every other character or factor column.
#' }
#'
#' An `add` column is appended to the result, taking the value `"old"` for
#' original rows and `"new"` for inserted padding rows.
#'
#' @param data    Data frame containing the trial layout.
#' @param pattern Character string of the form `"Var1:Var2"` naming the two
#'   spatial coordinate columns (rows and columns of the field grid).
#'   Default `"Row:Column"`.
#' @param match   Character vector of `Type` values that define the target
#'   sub-trial (e.g. test lines). The bounding box is computed from these
#'   plots only. Default `"DH"`.
#' @param split   Name of the blocking/grouping column. Extraction and padding
#'   are performed independently within each level. Default `"Block"`.
#' @param pad     Logical. If `TRUE` (default), missing grid cells within the
#'   bounding box are inserted as `NA` rows.
#' @param keep       Character vector of column names whose values should be
#'   carried into the padding rows (must be constant within each block).
#'   Defaults to `split`, which preserves the block identifier.
#' @param fill_value Character string written into every character or factor
#'   column of the padding rows, except for columns named in `keep` and the
#'   two spatial coordinate columns. Numeric columns always remain `NA`.
#'   Default `"Blank"`.
#' @param verbose Logical. If `TRUE`, prints a per-block message reporting the
#'   detected bounding box and the number of cells padded. Default `FALSE`.
#'
#' @return A data frame with the same columns as `data` plus an `add` column
#'   (`"old"` / `"new"`), ordered by the first then second spatial coordinate.
#'   Row and Column factors are re-levelled in ascending numeric order across
#'   the combined output.
#'
#' @examples
#' \dontrun{
#' # See scratch/demo_padTrial.R for a worked example
#' result <- padTrial(trial_df,
#'                    pattern    = "Row:Column",
#'                    match      = "DH",
#'                    split      = "Block",
#'                    fill_value = "Blank",
#'                    verbose    = TRUE)
#'
#' table(result$add)             # count original vs padded rows
#' subset(result, add == "new")  # inspect the inserted blank plots
#' }
#'
#' @export
padTrial <- function(data,
                     pattern    = "Row:Column",
                     match      = "DH",
                     split      = "Block",
                     pad        = TRUE,
                     keep       = split,
                     fill_value = "Blank",
                     verbose    = FALSE) {

  # ---- Input validation --------------------------------------------------
  pat <- unlist(strsplit(pattern, ":", fixed = TRUE))
  if (length(pat) != 2L)
    stop("'pattern' must be two column names separated by ':', e.g. \"Row:Column\".")

  if (!(split %in% names(data)))
    stop("'split' column \"", split, "\" not found in data.")

  if (!all(pat %in% names(data)))
    stop("Column(s) from 'pattern' not found in data: ",
         paste(setdiff(pat, names(data)), collapse = ", "), ".")

  if (!all(keep %in% names(data)))
    stop("'keep' column(s) not found in data: ",
         paste(setdiff(keep, names(data)), collapse = ", "), ".")

  if (!is.character(fill_value) || length(fill_value) != 1L)
    stop("'fill_value' must be a single character string.")

  if (!"Type" %in% names(data))
    stop("A 'Type' column is required to identify plot types via 'match'.")

  if (!any(as.character(data$Type) %in% match))
    stop("None of the 'match' value(s) found in the Type column: ",
         paste(match, collapse = ", "), ".")

  # ---- Internal helper ---------------------------------------------------
  # Safely coerce a factor or character column to numeric
  .as_num <- function(x) as.numeric(as.character(x))

  # ---- Process each block ------------------------------------------------
  block_list <- split(data, data[[split]])

  processed <- lapply(block_list, function(block) {

    block_id <- as.character(unique(block[[split]]))

    # Rows belonging to the target type
    target <- block[as.character(block$Type) %in% match, , drop = FALSE]

    if (nrow(target) == 0L) {
      if (verbose)
        message("Block ", block_id,
                ": no rows matching match = \"",
                paste(match, collapse = "/"), "\" -- block returned unchanged.")
      return(block)
    }

    # Bounding box
    row_range <- range(.as_num(target[[pat[1]]]))
    col_range <- range(.as_num(target[[pat[2]]]))

    if (verbose)
      message("Block ", block_id,
              ":  ", pat[1], " ", row_range[1], "-", row_range[2],
              "  |  ", pat[2], " ", col_range[1], "-", col_range[2])

    # Subset to bounding box
    in_bbox <- .as_num(block[[pat[1]]]) >= row_range[1] &
               .as_num(block[[pat[1]]]) <= row_range[2] &
               .as_num(block[[pat[2]]]) >= col_range[1] &
               .as_num(block[[pat[2]]]) <= col_range[2]

    sub <- droplevels(block[in_bbox, , drop = FALSE])
    sub$add <- "old"

    if (!pad) return(sub)

    # ---- Detect and fill missing cells -----------------------------------
    tabs    <- table(sub[[pat[1]]], sub[[pat[2]]])
    missing <- which(tabs == 0L, arr.ind = TRUE)

    if (nrow(missing) == 0L) {
      if (verbose) message("  -> no missing cells, no padding needed.")
      return(sub)
    }

    if (verbose)
      message("  -> padding ", nrow(missing), " missing cell(s).")

    # Build NA skeleton rows for padding
    n_pad <- nrow(missing)
    tp    <- sub[rep(1L, n_pad), , drop = FALSE]
    tp[]  <- lapply(tp, function(col) {
      if (is.factor(col))
        factor(rep(NA_character_, n_pad), levels = levels(col))
      else if (is.character(col))
        rep(NA_character_, n_pad)
      else
        rep(NA_real_, n_pad)
    })

    # Carry constant identifier columns into padding rows
    for (k in keep) {
      vals <- unique(sub[[k]])
      if (length(vals) != 1L)
        warning("'keep' column \"", k, "\" has ", length(vals),
                " distinct values in block ", block_id,
                " -- using the first value for padding rows.")
      tp[[k]] <- rep(vals[1L], n_pad)
    }

    # Set spatial coordinates of the padded cells
    tp[[pat[1]]] <- factor(rownames(tabs)[missing[, 1L]])
    tp[[pat[2]]] <- factor(colnames(tabs)[missing[, 2L]])

    # Fill all character/factor columns that are not keep or spatial coords
    fill_cols <- setdiff(
      names(tp)[vapply(tp, function(col) is.character(col) || is.factor(col),
                       logical(1L))],
      c(keep, pat, "add")
    )
    for (fc in fill_cols)
      tp[[fc]] <- if (is.factor(tp[[fc]])) {
        factor(rep(fill_value, n_pad),
               levels = union(levels(tp[[fc]]), fill_value))
      } else {
        rep(fill_value, n_pad)
      }

    tp$add      <- "new"
    rownames(tp) <- NULL

    rbind(sub, tp)
  })

  # ---- Combine blocks and re-order ---------------------------------------
  out          <- do.call(rbind, processed)
  rownames(out) <- NULL

  # Re-level spatial coordinates in ascending numeric order
  for (p in pat) {
    lvls    <- as.character(sort(unique(.as_num(out[[p]]))))
    out[[p]] <- factor(out[[p]], levels = lvls)
  }

  out[order(.as_num(out[[pat[1]]]), .as_num(out[[pat[2]]])), ]
}
