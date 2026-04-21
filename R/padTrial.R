#' Extract and Pad a Sub-Trial from a Field Trial Layout
#'
#' @description
#' Given a field trial data frame that mixes multiple plot types (e.g. test
#' lines, checks, guard rows), `padTrial()` extracts the rectangular sub-trial
#' occupied by a target plot type and, optionally, inserts placeholder rows for
#' any missing grid positions within that rectangle.
#'
#' **Steps performed per group (block):**
#' \enumerate{
#'   \item Check for duplicate Row × Column positions and stop with an
#'         informative error if any are found.
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
#' @param data     Data frame containing the trial layout.
#' @param pattern  Character string of the form `"Var1:Var2"` naming the two
#'   spatial coordinate columns (rows and columns of the field grid).
#'   Default `"Row:Column"`.
#' @param match    Character vector of values in `type_col` that define the
#'   target sub-trial (e.g. test lines). The bounding box is computed from
#'   these plots only. Default `"DH"`.
#' @param split    Character vector of one or more column names used to define
#'   independent groups (blocks) for extraction and padding. Pass `NULL` to
#'   treat the entire dataset as a single block. Default `"Block"`.
#' @param pad      Logical. If `TRUE` (default), missing grid cells within the
#'   bounding box are inserted as placeholder rows.
#' @param keep     Character vector of column names whose values should be
#'   carried into the padding rows (must be constant within each group).
#'   Defaults to `split`. Ignored when `split = NULL`.
#' @param fill_value Character string written into every character or factor
#'   column of the padding rows, except for columns named in `keep` and the
#'   two spatial coordinate columns. Numeric columns always remain `NA`.
#'   Default `"Blank"`.
#' @param type_col Name of the column in `data` that identifies plot types
#'   (e.g. `"DH"`, `"Check"`, `"Guard"`). Default `"Type"`.
#' @param verbose  Logical. If `TRUE`, prints a per-group message reporting the
#'   detected bounding box and the number of cells padded. Default `FALSE`.
#'
#' @return A data frame with the same columns as `data` plus an `add` column
#'   (`"old"` / `"new"`), ordered by the first then second spatial coordinate.
#'   The spatial coordinate columns are re-levelled as factors in ascending
#'   numeric order across the combined output. When multiple `split` columns
#'   are supplied the temporary grouping key column is dropped before returning.
#'
#' @examples
#' \dontrun{
#' # Single blocking column
#' result <- padTrial(trial_df,
#'                    pattern    = "Row:Column",
#'                    match      = "DH",
#'                    split      = "Block",
#'                    fill_value = "Blank",
#'                    verbose    = TRUE)
#'
#' # No blocking -- treat whole dataset as one group
#' result <- padTrial(trial_df, split = NULL)
#'
#' # Multi-environment: process each Site x Block independently
#' result <- padTrial(trial_df, split = c("Site", "Block"))
#'
#' # Custom type column name
#' result <- padTrial(trial_df, type_col = "PlotType", match = "Line")
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
                     type_col   = "Type",
                     verbose    = FALSE) {

  # ---- Input validation --------------------------------------------------
  pat <- unlist(strsplit(pattern, ":", fixed = TRUE))
  if (length(pat) != 2L)
    stop("'pattern' must be two column names separated by ':', ",
         "e.g. \"Row:Column\".")

  if (!all(pat %in% names(data)))
    stop("Column(s) from 'pattern' not found in data: ",
         paste(setdiff(pat, names(data)), collapse = ", "), ".")

  if (!is.null(split)) {
    if (!all(split %in% names(data)))
      stop("'split' column(s) not found in data: ",
           paste(setdiff(split, names(data)), collapse = ", "), ".")
    if (!all(keep %in% names(data)))
      stop("'keep' column(s) not found in data: ",
           paste(setdiff(keep, names(data)), collapse = ", "), ".")
  }

  if (!is.character(fill_value) || length(fill_value) != 1L)
    stop("'fill_value' must be a single character string.")

  if (!(type_col %in% names(data)))
    stop("'type_col' column \"", type_col, "\" not found in data.")

  if (!any(as.character(data[[type_col]]) %in% match))
    stop("None of the 'match' value(s) found in column \"", type_col, "\": ",
         paste(match, collapse = ", "), ".")

  # ---- Internal helper ---------------------------------------------------
  .as_num <- function(x) as.numeric(as.character(x))

  # ---- Build grouping key ------------------------------------------------
  # split = NULL  -> single synthetic group
  # split = one column -> use directly
  # split = multiple columns -> paste into a temporary composite key
  GRP_COL <- ".padTrial_grp"

  if (is.null(split)) {
    data[[GRP_COL]] <- "all"
    grp_keep        <- character(0L)
  } else if (length(split) == 1L) {
    data[[GRP_COL]] <- as.character(data[[split]])
    grp_keep        <- keep
  } else {
    data[[GRP_COL]] <- do.call(paste, c(lapply(split, function(s)
      as.character(data[[s]])), list(sep = "\u00b7")))   # middle-dot separator
    grp_keep        <- keep
  }

  # ---- Process each group ------------------------------------------------
  group_list <- split(data, data[[GRP_COL]])

  processed <- lapply(group_list, function(block) {

    grp_id <- as.character(unique(block[[GRP_COL]]))
    lbl    <- if (is.null(split)) "All data" else paste("Group", grp_id)

    # -- Duplicate Row x Column check ------------------------------------
    coord_key <- paste(.as_num(block[[pat[1]]]),
                       .as_num(block[[pat[2]]]), sep = ",")
    dups <- coord_key[duplicated(coord_key)]
    if (length(dups) > 0L) {
      dup_str <- paste(unique(dups), collapse = "; ")
      stop(lbl, ": duplicate Row \u00d7 Column position(s) found: ",
           dup_str, ".\n",
           "  Each grid cell must appear at most once per group. ",
           "Check your data for repeated records.")
    }

    # -- Target plots and bounding box -----------------------------------
    target <- block[as.character(block[[type_col]]) %in% match, , drop = FALSE]

    if (nrow(target) == 0L) {
      if (verbose)
        message(lbl, ": no rows matching match = \"",
                paste(match, collapse = "/"), "\" -- group returned unchanged.")
      block[[GRP_COL]] <- NULL
      return(block)
    }

    row_range <- range(.as_num(target[[pat[1]]]))
    col_range <- range(.as_num(target[[pat[2]]]))

    if (verbose)
      message(lbl, ":  ", pat[1], " ", row_range[1], "-", row_range[2],
              "  |  ", pat[2], " ", col_range[1], "-", col_range[2])

    # -- Subset to bounding box ------------------------------------------
    in_bbox <- .as_num(block[[pat[1]]]) >= row_range[1] &
               .as_num(block[[pat[1]]]) <= row_range[2] &
               .as_num(block[[pat[2]]]) >= col_range[1] &
               .as_num(block[[pat[2]]]) <= col_range[2]

    sub <- droplevels(block[in_bbox, , drop = FALSE])
    sub[[GRP_COL]] <- NULL   # drop temporary grouping key
    sub$add <- "old"

    if (!pad) return(sub)

    # -- Detect and fill missing cells -----------------------------------
    tabs    <- table(sub[[pat[1]]], sub[[pat[2]]])
    missing <- which(tabs == 0L, arr.ind = TRUE)

    if (nrow(missing) == 0L) {
      if (verbose) message("  -> no missing cells, no padding needed.")
      return(sub)
    }

    if (verbose)
      message("  -> padding ", nrow(missing), " missing cell(s).")

    # Build NA skeleton rows
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
    for (k in grp_keep) {
      vals <- unique(sub[[k]])
      if (length(vals) != 1L)
        warning("'keep' column \"", k, "\" has ", length(vals),
                " distinct values in ", lbl,
                " -- using the first value for padding rows.")
      tp[[k]] <- rep(vals[1L], n_pad)
    }

    # Set spatial coordinates of the padded cells
    tp[[pat[1]]] <- factor(rownames(tabs)[missing[, 1L]])
    tp[[pat[2]]] <- factor(colnames(tabs)[missing[, 2L]])

    # Fill all remaining character/factor columns with fill_value
    fill_cols <- setdiff(
      names(tp)[vapply(tp, function(col) is.character(col) || is.factor(col),
                       logical(1L))],
      c(grp_keep, pat, "add")
    )
    for (fc in fill_cols)
      tp[[fc]] <- if (is.factor(tp[[fc]])) {
        factor(rep(fill_value, n_pad),
               levels = union(levels(tp[[fc]]), fill_value))
      } else {
        rep(fill_value, n_pad)
      }

    tp$add       <- "new"
    rownames(tp) <- NULL

    rbind(sub, tp)
  })

  # ---- Combine groups and re-order ---------------------------------------
  out          <- do.call(rbind, processed)
  rownames(out) <- NULL

  # Re-level spatial coordinates in ascending numeric order
  for (p in pat) {
    lvls    <- as.character(sort(unique(.as_num(out[[p]]))))
    out[[p]] <- factor(out[[p]], levels = lvls)
  }

  out[order(.as_num(out[[pat[1]]]), .as_num(out[[pat[2]]])), ]
}
