#' Multiple Comparison Criteria for ASReml-R Predicted Values
#'
#' @description
#' Computes a pairwise comparison criterion (HSD, LSD, or Bonferroni-corrected
#' LSD) for predicted values obtained from an ASReml-R model, optionally
#' within groups defined by a subset of the factors in `term`.
#'
#' For each group the function computes the **average SED** across all
#' \eqn{\binom{n_g}{2}} pairs of predictions:
#'
#' \deqn{
#'   \overline{\text{SED}} = \frac{1}{\binom{n_g}{2}}
#'   \sum_{i < j} \sqrt{V_{ii} + V_{jj} - 2V_{ij}}
#' }
#'
#' where \eqn{V_{ij}} is the \eqn{(i,j)} element of the prediction error
#' variance-covariance matrix.  The comparison criterion is then:
#'
#' \describe{
#'   \item{`"HSD"` (Tukey)}{
#'     \eqn{\displaystyle\frac{\overline{\text{SED}}}{\sqrt{2}} \times
#'     q_{1-\alpha}(n_g,\, \nu)}, where \eqn{q} is the Studentised range
#'     quantile with \eqn{n_g} means and \eqn{\nu = \texttt{model\$nedf}}
#'     error df.  Two means are declared significantly different if their
#'     absolute difference exceeds the HSD.}
#'   \item{`"LSD"` (unadjusted)}{
#'     \eqn{\overline{\text{SED}} \times t_{\alpha/2,\,\nu}}.
#'     No multiplicity correction; liberal for large numbers of comparisons.}
#'   \item{`"Bonferroni"`}{
#'     \eqn{\overline{\text{SED}} \times
#'     t_{\alpha/(2m),\,\nu}}, where \eqn{m = \binom{n_g}{2}} is the
#'     number of pairs.  Controls the family-wise error rate.}
#' }
#'
#' @param model   An ASReml-R V4 model object.
#' @param term    Character string specifying the `classify` term, e.g.
#'   `"Treatment:Genotype"`.  Predictions are obtained via
#'   `asreml::predict.asreml()`.
#' @param by      Grouping specification.  Comparisons are performed
#'   independently within each level of the group.  Accepted forms:
#'   \describe{
#'     \item{`NULL` (default)}{All predictions form one group.}
#'     \item{Single string `"Site"`}{Split by the `Site` factor in `term`.}
#'     \item{Colon string `"Site:Treatment"`}{Combine `Site` and `Treatment`
#'       into a single grouping variable and split by their interaction.}
#'     \item{Character vector `c("Site", "Treatment")`}{Equivalent to the
#'       colon form above.}
#'   }
#'   All `by` variables must be present in `term`.  At least one variable in
#'   `term` must remain for comparisons to be made.
#' @param type    Comparison type.  One of `"HSD"` (Tukey, default),
#'   `"LSD"` (unadjusted), or `"Bonferroni"`.
#' @param pev     Logical.  If `TRUE` (default) the prediction error
#'   variance (PEV) from `predict()$vcov` is used.  If `FALSE` and `term`
#'   is in the random part of the model, the posterior variance
#'   \eqn{\bm{G} - \text{PEV}} is used instead.
#' @param alpha   Significance level.  Default `0.05`.
#' @param ...     Additional arguments forwarded to
#'   `asreml::predict.asreml()`.
#'
#' @return A data frame with one row per non-missing predicted value,
#'   containing:
#'   \describe{
#'     \item{Factor columns}{All factor variables from `term`.}
#'     \item{`predicted.value`}{The BLUE or BLUP.}
#'     \item{`std.error`}{Prediction standard error.}
#'     \item{`{type}`}{The comparison criterion (constant within each group).
#'       Two predicted values in the same group are significantly different if
#'       their absolute difference exceeds this value.}
#'     \item{`avsed`}{The average SED for the group (constant within group).}
#'   }
#'
#' @examples
#' \dontrun{
#' ## HSD for Treatment:Genotype, one group per Site
#' res <- compare(model,
#'                term = "Treatment:Site:Genotype",
#'                by   = "Site",
#'                type = "HSD")
#'
#' ## LSD computed within each Treatment x Site combination
#' res2 <- compare(model,
#'                 term = "Treatment:Site:Genotype",
#'                 by   = c("Treatment", "Site"),
#'                 type = "LSD")
#'
#' ## Equivalently using colon notation
#' res3 <- compare(model,
#'                 term = "Treatment:Site:Genotype",
#'                 by   = "Treatment:Site",
#'                 type = "LSD")
#'
#' ## Bonferroni-corrected LSD, single global group
#' res4 <- compare(model,
#'                 term  = "Treatment:Genotype",
#'                 type  = "Bonferroni",
#'                 alpha = 0.05)
#'
#' ## Posterior variance (BLUPs)
#' res5 <- compare(model,
#'                 term = "Treatment:Genotype",
#'                 by   = "Treatment",
#'                 pev  = FALSE)
#' }
#'
#' @seealso [waldTest()] for pairwise p-values, `asreml::predict.asreml()`
#' @export
compare <- function(model, term, by = NULL,
                    type  = c("HSD", "LSD", "Bonferroni"),
                    pev   = TRUE,
                    alpha = 0.05,
                    ...) {

  type <- match.arg(type)

  if (!inherits(model, "asreml"))
    stop("'model' must be an object of class \"asreml\".")
  if (alpha <= 0 || alpha >= 1)
    stop("'alpha' must be strictly between 0 and 1.")

  # ---- Predictions and variance matrix ----------------------------------
  pred  <- predict(model, classify = term, vcov = TRUE, ...)
  terms <- unlist(strsplit(term, ":"))
  pv    <- pred$pvals

  # Remove NA predictions; align vcov
  whna  <- !is.na(pv$predicted.value)
  pv    <- pv[whna, , drop = FALSE]
  Vcov  <- as.matrix(pred$vcov)[whna, whna, drop = FALSE]
  rownames(pv) <- NULL

  # Posterior variance (G - PEV) when pev = FALSE
  if (!pev) {
    rterms <- all.vars(model$call$random)
    if (!all(terms %in% rterms)) {
      warning("'pev = FALSE' requires all terms to be in the random formula. ",
              "Falling back to PEV.", call. = FALSE)
    } else {
      varm <- .asreml_vparams(model, term)
      if (is.null(varm))
        warning("Could not retrieve G-matrix for '", term,
                "'. Falling back to PEV.", call. = FALSE)
      else {
        len  <- if (length(terms) > 1L)
                  as.integer(table(pv[, 1L])[[1L]])
                else
                  nrow(pv)
        Vcov <- (kronecker(varm, diag(len)) - as.matrix(pred$vcov))[whna, whna,
                                                                     drop = FALSE]
      }
    }
  }

  # ---- Parse 'by' grouping ----------------------------------------------
  by_vars <- NULL
  if (!is.null(by)) {

    # Accept: "Site", "Site:Treatment", or c("Site","Treatment")
    by_vars <- if (length(by) == 1L)
      unlist(strsplit(by, ":", fixed = TRUE))
    else
      as.character(by)

    bad <- setdiff(by_vars, terms)
    if (length(bad))
      stop("'by' variables not found in 'term': ",
           paste(bad, collapse = ", "), ".")

    remaining <- setdiff(terms, by_vars)
    if (!length(remaining))
      stop("'by' accounts for all variables in 'term'  -- ",
           "no levels remain within groups to compare.")

    # Construct combined group label (single string per row)
    pv[["grp__"]] <- if (length(by_vars) == 1L)
      as.character(pv[[by_vars]])
    else
      apply(pv[, by_vars, drop = FALSE], 1L, paste, collapse = ":")

  } else {
    pv[["grp__"]] <- "All"
  }

  um     <- unique(pv[["grp__"]])
  df_err <- model$nedf

  if (is.null(df_err) || is.na(df_err))
    stop("'model$nedf' is NULL or NA  -- cannot compute comparison criterion.")

  # ---- Compute criterion per group --------------------------------------
  crit_vec <- numeric(nrow(pv))
  sed_vec  <- numeric(nrow(pv))

  for (i in seq_along(um)) {

    gind <- which(pv[["grp__"]] == um[i])
    n_g  <- length(gind)

    if (n_g < 2L) {
      warning("Group '", um[i], "' has only one observation  -- skipped.",
              call. = FALSE)
      crit_vec[gind] <- NA_real_
      sed_vec[gind]  <- NA_real_
      next
    }

    svar <- Vcov[gind, gind, drop = FALSE]

    # All C(n_g, 2) pairwise SED^2 values
    dv   <- diag(svar)
    sed2 <- outer(dv, dv, "+") - 2 * svar     # n_g x n_g matrix of SED^2
    sed2 <- sed2[lower.tri(sed2)]               # lower triangle only
    sed2[sed2 < 0] <- NA_real_                  # guard against rounding below 0

    # Arithmetic mean of SEDs (not RMS)
    seds    <- sqrt(sed2)
    avsed   <- mean(seds, na.rm = TRUE)
    m_pairs <- sum(!is.na(seds))

    # Comparison criterion
    crit <- switch(type,
      HSD        = (avsed / sqrt(2)) *
                     qtukey(1 - alpha, nmeans = n_g, df = df_err),
      LSD        =  avsed *
                     qt(alpha / 2, df = df_err, lower.tail = FALSE),
      Bonferroni =  avsed *
                     qt(alpha / (2 * m_pairs), df = df_err, lower.tail = FALSE)
    )

    crit_vec[gind] <- crit
    sed_vec[gind]  <- avsed
  }

  # ---- Assemble output --------------------------------------------------
  # Keep all factor columns from 'term' plus predicted.value and std.error
  keep <- intersect(c(terms, "predicted.value", "std.error"), names(pv))
  out  <- pv[, keep, drop = FALSE]

  out[[type]]    <- crit_vec
  out[["avsed"]] <- sed_vec
  rownames(out)  <- NULL
  out
}
