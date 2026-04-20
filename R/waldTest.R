# ---- Private helpers ----------------------------------------------------

## Construct row labels from the factor columns of a pvals data frame.
## Excludes the 'status' column and optionally one grouping variable.
#' @noRd
.rowLabels <- function(pvals, exclude = NULL) {
  fac_cols <- names(pvals)[vapply(pvals, is.factor, logical(1L))]
  fac_cols <- setdiff(fac_cols, c("status", exclude))
  if (!length(fac_cols))
    return(as.character(seq_len(nrow(pvals))))
  if (length(fac_cols) == 1L)
    return(as.character(pvals[[fac_cols]]))
  apply(pvals[, fac_cols, drop = FALSE], 1L, paste, collapse = ":")
}

## Build the C(k,2) pairwise contrast matrix for k levels.
#' @noRd
.pairwiseMat <- function(k) {
  pairs <- combn(k, 2L)
  mat   <- matrix(0, ncol(pairs), k)
  for (i in seq_len(ncol(pairs))) {
    mat[i, pairs[1L, i]] <-  1
    mat[i, pairs[2L, i]] <- -1
  }
  mat
}

## Resolve coef specification to integer positions within a group.
## 'coef'   : integer indices or character labels
## 'labels' : character vector of row labels for the current group
#' @noRd
.resolveCoef <- function(coef, labels) {
  if (is.character(coef)) {
    idx <- match(coef, labels)
    bad <- coef[is.na(idx)]
    if (length(bad))
      stop("Labels not found in predictions: ", paste(bad, collapse = ", "))
    idx
  } else {
    if (any(coef < 1L) || any(coef > length(labels)))
      stop("Coefficient index out of bounds (group has ", length(labels),
           " rows).")
    as.integer(coef)
  }
}

## Build a contrast matrix row (or block of rows) for one 'con' element.
## Returns a list: mat (nrow x n_g), row_labels (nrow x 2 character matrix).
#' @noRd
.buildCRows <- function(el, labels) {
  n_g    <- length(labels)
  coef   <- .resolveCoef(el$coef, labels)
  k      <- length(coef)

  ## Resolve comp to a matrix
  if (is.null(el$comp) || identical(el$comp, "helmert")) {
    if (k < 2L) stop("Need at least 2 coefficients for a contrast.")
    comp <- matrix(contr.helmert(k)[, k - 1L], nrow = 1L)
    warning("Default Helmert contrast used for: ",
            paste(labels[coef], collapse = ", "),
            call. = FALSE)
  } else if (identical(el$comp, "pairwise")) {
    comp <- .pairwiseMat(k)
  } else if (is.matrix(el$comp)) {
    if (ncol(el$comp) != k)
      stop("Contrast matrix must have ", k, " columns (length of 'coef').")
    comp <- el$comp
  } else {
    if (length(el$comp) != k)
      stop("Contrast vector must have length ", k, " (length of 'coef').")
    comp <- matrix(el$comp, nrow = 1L)
  }

  nr  <- nrow(comp)
  mat <- matrix(0, nr, n_g)
  mat[, coef] <- comp

  ## Auto-generate row names (left vs right based on sign of coefficients)
  rnam <- matrix("", nr, 2L)
  for (i in seq_len(nr)) {
    rnam[i, 1L] <- paste(labels[coef[comp[i, ] > 0]], collapse = ":")
    rnam[i, 2L] <- paste(labels[coef[comp[i, ] < 0]], collapse = ":")
  }

  ## Override with user-supplied group names
  if (!is.null(el$group)) {
    grp <- el$group
    if (!all(names(grp) %in% c("left", "right")))
      stop("'group' must have names 'left' and/or 'right'.")
    if (!is.null(grp$left)) {
      gL <- rep_len(grp$left, nr)
      rnam[!is.na(gL), 1L] <- gL[!is.na(gL)]
    }
    if (!is.null(grp$right)) {
      gR <- rep_len(grp$right, nr)
      rnam[!is.na(gR), 2L] <- gR[!is.na(gR)]
    }
  }

  list(mat = mat, rnam = rnam)
}

## Build a selection matrix for one 'zero' element.
#' @noRd
.buildZRows <- function(el, labels) {
  n_g  <- length(labels)
  coef <- .resolveCoef(el$coef, labels)
  q    <- length(coef)
  Z    <- matrix(0, q, n_g)
  for (j in seq_len(q)) Z[j, coef[j]] <- 1L

  test_name <- if (!is.null(el$group)) el$group else
    paste(labels[coef], collapse = ":")

  list(Z = Z, name = test_name, df = q)
}


# ---- Main function (works on predict output) ----------------------------

#' Wald Tests for Fixed-Effect Contrasts Using Predicted Values
#'
#' @description
#' Performs Wald tests on linear contrasts of predicted values obtained from
#' [asreml::predict.asreml()].  Working directly with `predict()` output
#' means contrasts are specified through **meaningful factor-level labels**
#' rather than raw coefficient indices, and the prediction error
#' variance--covariance matrix (`pred$vcov`) is used directly  -- no access to
#' model internals required.
#'
#' Two test types are available for each element of `cc`:
#' \describe{
#'   \item{`"con"`  -- contrasts}{Test \eqn{H_0: \bm{c}^\top\hat{\bm{\tau}} = 0}
#'     for a specified contrast vector \eqn{\bm{c}}.  The Wald statistic is
#'     \eqn{W = (\bm{c}^\top\hat{\bm{\tau}})^2 / \bm{c}^\top\text{Vcov}\,\bm{c}
#'     \sim \chi^2_1} (or \eqn{F_{1,\nu}} with `test = "F"`).  Set
#'     `comp = "pairwise"` to automatically generate all \eqn{\binom{k}{2}}
#'     pairwise contrasts for \eqn{k} specified levels.}
#'   \item{`"zero"`  -- joint zero equality}{Test \eqn{H_0: \bm{Z}\hat{\bm{\tau}} = \bm{0}}
#'     for a set of \eqn{q} coefficients simultaneously.  The statistic is
#'     \eqn{W = \hat{\bm{\tau}}^\top\bm{Z}^\top(\bm{Z}\,\text{Vcov}\,\bm{Z}^\top)^{-1}
#'     \bm{Z}\hat{\bm{\tau}} \sim \chi^2_q} (or \eqn{F_{q,\nu}/q}).}
#' }
#'
#' @param pred     List returned by `predict(model, classify = ..., vcov = TRUE)`.
#'   Must contain elements `pvals` (data frame of predictions) and `vcov`
#'   (prediction error variance--covariance matrix).
#' @param cc       Named list of test specifications.  Each element is itself
#'   a list with fields:
#'   \describe{
#'     \item{`coef`}{Integer indices **or** character labels identifying which
#'       rows of `pred$pvals` are involved.  When `by` is used, indices/labels
#'       refer to rows **within the group**.  Character labels are matched
#'       against the pasted factor-level combination (e.g. `"T0:G01"` for a
#'       `Treatment:Genotype` prediction).}
#'     \item{`type`}{`"con"` (contrast) or `"zero"` (joint zero test).}
#'     \item{`comp`}{For `type = "con"` only. One of:
#'       a numeric contrast vector (length = `length(coef)`);
#'       a numeric contrast matrix (rows = contrasts, cols = `length(coef)`);
#'       `"pairwise"` to generate all pairwise differences automatically;
#'       `NULL` or `"helmert"` for the default last Helmert contrast
#'       (a warning is issued).}
#'     \item{`group`}{Optional named list `list(left = ..., right = ...)` to
#'       override the auto-generated `"A vs B"` row labels.}
#'   }
#' @param by       Character string **or character vector** naming one or more
#'   factor columns in `pred$pvals` to split by.  When a single name is
#'   supplied the rows are split by that column's levels.  When a vector of
#'   names is supplied (e.g. `by = c("Site", "Year")`) the levels of those
#'   columns are pasted together to form an interaction grouping variable, and
#'   the output grouping column is named with the same pasted string (e.g.
#'   `"Site:Year"`).  The same `cc` specification is applied independently
#'   within each group, with `coef` referring to within-group row positions.
#'   Default `NULL` (no splitting).
#' @param test     `"Wald"` (default) for asymptotic \eqn{\chi^2} tests, or
#'   `"F"` for \eqn{F}-tests.  `"F"` requires `df_error`.
#' @param df_error Numeric error degrees of freedom for F-tests (e.g.
#'   `model$nedf`).  Ignored when `test = "Wald"`.
#' @param adjust   P-value adjustment method passed to [p.adjust()].  Applied
#'   across all contrasts within each group.  Options:
#'   \describe{
#'     \item{`"none"`}{No adjustment (default).}
#'     \item{`"bonferroni"`}{Bonferroni correction  -- controls family-wise error
#'       rate (FWER) conservatively.}
#'     \item{`"holm"`}{Holm step-down  -- controls FWER less conservatively than
#'       Bonferroni.}
#'     \item{`"fdr"` / `"BH"`}{Benjamini--Hochberg False Discovery Rate.
#'       Controls the expected proportion of false discoveries among rejected
#'       hypotheses.  Best choice when many contrasts are tested and some true
#'       effects are expected.}
#'     \item{`"BY"`}{Benjamini--Yekutieli FDR.  Valid under arbitrary (including
#'       positive) dependence between test statistics; more conservative than
#'       `"BH"` but provides stronger guarantees.}
#'   }
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`Contrasts`}{Data frame of contrast test results (or `NULL` if no
#'     `"con"` elements were supplied).  Columns: grouping variable (if `by`
#'     used), `Comparison`, `Estimate`, `Std.Error`,
#'     `Wald.Statistic` or `F.Statistic`, `df` (numerator), `P.Value`.}
#'   \item{`Zero`}{Data frame of joint zero test results (or `NULL`).
#'     Columns: grouping variable (if `by` used), `Test`,
#'     `Wald.Statistic` or `F.Statistic`, `df`, `P.Value`.}
#'   \item{`test`}{The `test` argument used.}
#'   \item{`adjust`}{The `adjust` argument used.}
#' }
#'
#' @examples
#' \dontrun{
#' ## --- predict first, then test ---
#' pred <- predict(model, classify = "Treatment", vcov = TRUE)
#'
#' ## Pairwise contrasts among all three treatment levels
#' res <- waldTest(pred,
#'                 cc = list(list(coef = c("N0","N1","N2"),
#'                                type = "con",
#'                                comp = "pairwise")))
#'
#' ## Manual contrast: N2 vs N0
#' res2 <- waldTest(pred,
#'                  cc = list(list(coef = c("N0","N2"),
#'                                 type = "con",
#'                                 comp = c(-1, 1),
#'                                 group = list(left="N2", right="N0"))))
#'
#' ## Joint zero test: are N1 and N2 effects simultaneously zero?
#' res3 <- waldTest(pred,
#'                  cc = list(list(coef = c("N1","N2"),
#'                                 type = "zero",
#'                                 group = "N1:N2 joint")))
#'
#' ## F-tests with denominator df from model
#' res4 <- waldTest(pred, cc = list(...),
#'                  test = "F", df_error = model$nedf)
#'
#' ## Run within each site (by-group, single factor)
#' pred_site <- predict(model, classify = "Treatment:Site", vcov = TRUE)
#' res5 <- waldTest(pred_site,
#'                  cc   = list(list(coef = c("N0","N1","N2"),
#'                                   type = "con", comp = "pairwise")),
#'                  by   = "Site",
#'                  test = "F", df_error = model$nedf,
#'                  adjust = "fdr")
#'
#' ## Run within each Site x Year combination (multi-factor by)
#' pred_sy <- predict(model, classify = "Treatment:Site:Year", vcov = TRUE)
#' res5b <- waldTest(pred_sy,
#'                   cc  = list(list(coef = c("N0","N1","N2"),
#'                                   type = "con", comp = "pairwise")),
#'                   by  = c("Site", "Year"),
#'                   test = "F", df_error = model$nedf)
#'
#' ## Convenience: pass the model directly
#' res6 <- waldTest(model, classify = "Treatment",
#'                  cc = list(list(coef=c("N0","N2"), type="con", comp=c(-1,1))),
#'                  test = "F")
#' }
#'
#' @seealso [asreml::predict.asreml()]
#' @export
waldTest <- function(pred, cc, by = NULL,
                     test     = c("Wald", "F"),
                     df_error = NULL,
                     adjust   = c("none","bonferroni","holm","fdr","BH","BY")) {

  test   <- match.arg(test)
  adjust <- match.arg(adjust)

  # ---- Validate pred input -----------------------------------------------
  if (!is.list(pred) || !all(c("pvals","vcov") %in% names(pred)))
    stop("'pred' must be the list returned by predict(model, vcov = TRUE).")

  pvals <- pred$pvals
  Vcov  <- as.matrix(pred$vcov)

  # Remove NA predictions and align vcov
  whna  <- !is.na(pvals$predicted.value)
  pvals <- pvals[whna, , drop = FALSE]
  Vcov  <- Vcov[whna, whna, drop = FALSE]
  rownames(pvals) <- NULL

  if (nrow(pvals) == 0L)
    stop("No non-missing predicted values in 'pred'.")

  if (test == "F" && is.null(df_error))
    stop("'df_error' must be supplied when test = \"F\".")

  # ---- Validate cc -------------------------------------------------------
  for (i in seq_along(cc)) {
    el <- cc[[i]]
    if (!all(names(el) %in% c("coef","type","comp","group")))
      stop("cc[[", i, "]] has unrecognised fields.  ",
           "Allowed: 'coef', 'type', 'comp', 'group'.")
    if (!el$type %in% c("con","zero"))
      stop("cc[[", i, "]]$type must be \"con\" or \"zero\".")
    if (el$type == "con" && is.null(el$comp))
      warning("cc[[", i, "]]$comp is NULL  -- default Helmert contrast will be used.",
              call. = FALSE)
  }

  ctype <- vapply(cc, `[[`, character(1L), "type")
  cons  <- cc[ctype == "con"]
  zeros <- cc[ctype == "zero"]

  # ---- Set up grouping ---------------------------------------------------
  if (!is.null(by)) {
    missing_by <- setdiff(by, names(pvals))
    if (length(missing_by))
      stop("'by' variable(s) not found in predictions: ",
           paste(missing_by, collapse = ", "))
    if (length(by) == 1L) {
      by_vals    <- as.character(pvals[[by]])
      by_label   <- by          # column name used in output
      by_exclude <- by          # column(s) excluded from .rowLabels
    } else {
      by_vals    <- apply(pvals[, by, drop = FALSE], 1L, paste, collapse = ":")
      by_label   <- paste(by, collapse = ":")
      by_exclude <- by
    }
  } else {
    by_vals    <- rep("All", nrow(pvals))
    by_label   <- "Group"
    by_exclude <- "Group"
  }
  um <- unique(by_vals)

  # ---- Main loop: one iteration per group --------------------------------
  con_list  <- vector("list", length(um))
  zero_list <- vector("list", length(um))

  for (gi in seq_along(um)) {

    g_idx  <- which(by_vals == um[gi])
    tau_g  <- pvals$predicted.value[g_idx]
    Vcov_g <- Vcov[g_idx, g_idx, drop = FALSE]
    n_g    <- length(g_idx)

    # Row labels within this group (exclude the 'by' variable(s))
      labels_g <- .rowLabels(pvals[g_idx, , drop = FALSE], exclude = by_exclude)

    # ---- Contrast tests --------------------------------------------------
    if (length(cons)) {

      Crow_list  <- lapply(cons, .buildCRows, labels = labels_g)
      Cmat_full  <- do.call(rbind, lapply(Crow_list, `[[`, "mat"))
      Cnam_full  <- do.call(rbind, lapply(Crow_list, `[[`, "rnam"))
      nr         <- nrow(Cmat_full)

      wald_c <- se_c <- est_c <- numeric(nr)
      for (i in seq_len(nr)) {
        c_i       <- Cmat_full[i, ]
        cv        <- drop(c_i %*% Vcov_g %*% c_i)
        se_c[i]   <- sqrt(cv)
        est_c[i]  <- drop(c_i %*% tau_g)
        wald_c[i] <- est_c[i]^2 / cv
      }

      pval_c <- if (test == "F")
        pf(wald_c, 1L, df_error, lower.tail = FALSE)
      else
        pchisq(wald_c, 1L, lower.tail = FALSE)

      if (adjust != "none") pval_c <- p.adjust(pval_c, method = adjust)

      stat_col <- if (test == "F") "F.Statistic" else "Wald.Statistic"
      cdf <- data.frame(
        .grp       = um[gi],
        Comparison = paste(Cnam_full[, 1L], Cnam_full[, 2L], sep = " vs "),
        Estimate   = round(est_c,  6L),
        Std.Error  = round(se_c,   6L),
        round(wald_c, 6L),
        df         = 1L,
        P.Value    = round(pval_c, 6L),
        stringsAsFactors = FALSE
      )
      names(cdf)[names(cdf) == ".grp"]              <- by_label
      names(cdf)[names(cdf) == "round.wald_c..6L."] <- stat_col
      con_list[[gi]] <- cdf
    }

    # ---- Zero equality tests ---------------------------------------------
    if (length(zeros)) {

      zrow_list <- lapply(zeros, .buildZRows, labels = labels_g)
      nr_z      <- length(zrow_list)
      wald_z    <- df_z <- numeric(nr_z)
      znam      <- character(nr_z)
      pval_z    <- numeric(nr_z)

      for (j in seq_len(nr_z)) {
        Z    <- zrow_list[[j]]$Z
        q    <- zrow_list[[j]]$df
        ctau <- drop(Z %*% tau_g)
        varZ <- Z %*% Vcov_g %*% t(Z)

        wstat <- tryCatch(
          drop(ctau %*% solve(varZ) %*% ctau),
          error = function(e) {
            warning("Singular variance matrix in zero test '",
                    zrow_list[[j]]$name, "'. Returning NA.", call. = FALSE)
            NA_real_
          }
        )

        wald_z[j] <- if (test == "F") wstat / q else wstat
        df_z[j]   <- q
        pval_z[j] <- if (is.na(wstat)) NA_real_ else
          if (test == "F") pf(wald_z[j], q, df_error, lower.tail = FALSE)
          else             pchisq(wstat,  q, lower.tail = FALSE)
        znam[j]   <- zrow_list[[j]]$name
      }

      if (adjust != "none") pval_z <- p.adjust(pval_z, method = adjust)

      stat_col <- if (test == "F") "F.Statistic" else "Wald.Statistic"
      zdf <- data.frame(
        .grp    = um[gi],
        Test    = znam,
        round(wald_z, 6L),
        df      = df_z,
        P.Value = round(pval_z, 6L),
        stringsAsFactors = FALSE
      )
      names(zdf)[names(zdf) == ".grp"]              <- by_label
      names(zdf)[names(zdf) == "round.wald_z..6L."] <- stat_col
      zero_list[[gi]] <- zdf
    }
  }

  # ---- Assemble and return -----------------------------------------------
  con_out  <- do.call(rbind, con_list[!vapply(con_list,  is.null, logical(1L))])
  zero_out <- do.call(rbind, zero_list[!vapply(zero_list, is.null, logical(1L))])

  # Drop the 'Group' column if no by was originally supplied
  if (identical(by_label, "Group")) {
    if (!is.null(con_out))  con_out$Group  <- NULL
    if (!is.null(zero_out)) zero_out$Group <- NULL
  }

  rownames(con_out)  <- NULL
  rownames(zero_out) <- NULL

  invisible(list(
    Contrasts = if (!is.null(con_out)  && nrow(con_out)  > 0L) con_out  else NULL,
    Zero      = if (!is.null(zero_out) && nrow(zero_out) > 0L) zero_out else NULL,
    test      = test,
    adjust    = adjust
  ))
}


# ---- Convenience S3 method for asreml models ----------------------------

#' @describeIn waldTest Convenience method that calls [asreml::predict.asreml()]
#'   internally.  `classify` is required; `test = "F"` uses `object$nedf`
#'   automatically.
#'
#' @param object   An ASReml-R V4 model object.
#' @param classify Character string passed to [asreml::predict.asreml()] as
#'   the `classify` argument.
#' @param ...      Additional arguments forwarded to
#'   [asreml::predict.asreml()].
#'
#' @export
waldTest.asreml <- function(object, classify, cc, by = NULL,
                             test     = c("Wald", "F"),
                             adjust   = c("none","bonferroni","holm","fdr","BH","BY"),
                             ...) {

  if (!inherits(object, "asreml"))
    stop("'object' must be of class \"asreml\".")

  test   <- match.arg(test)
  adjust <- match.arg(adjust)

  pred <- predict(object, classify = classify, vcov = TRUE, ...)

  df_error <- if (test == "F") {
    if (is.null(object$nedf))
      stop("'object$nedf' is NULL  -- cannot use F-test. Use test = \"Wald\".")
    object$nedf
  } else {
    NULL
  }

  waldTest(pred       = pred,
           cc         = cc,
           by         = by,
           test       = test,
           df_error   = df_error,
           adjust     = adjust)
}


# ---- print method -------------------------------------------------------

#' @export
print.waldTest <- function(x, digits = 4L, ...) {
  if (!is.null(x$Contrasts)) {
    cat("\nContrast Tests (", x$test, "):\n", sep = "")
    print(format(x$Contrasts, digits = digits), row.names = FALSE)
  }
  if (!is.null(x$Zero)) {
    cat("\nJoint Zero-Equality Tests (", x$test, "):\n", sep = "")
    print(format(x$Zero, digits = digits), row.names = FALSE)
  }
  if (!is.null(x$adjust) && x$adjust != "none")
    cat("\nP-value adjustment: ", x$adjust, "\n", sep = "")
  invisible(x)
}
