#' Multivariate Fixed-Effects Regression of Treatment BLUEs Within Groups
#'
#' @description
#' The fixed-effects (BLUE) analogue of [randomRegress()].  For each group
#' defined by the `by` argument, the function regresses the BLUEs of each
#' conditioned treatment on the BLUEs of its conditioning set using ordinary
#' least squares, returning a **response index** (regression residual) for
#' every genotype.
#'
#' For treatment \eqn{j} with conditioning set \eqn{A_j} of size
#' \eqn{a = |A_j|}, the OLS fit is:
#'
#' \deqn{
#'   \hat{\bm{\tau}}_j = \bm{X}_j \hat{\bm{\beta}}_j + \tilde{\bm{\tau}}_j,
#'   \qquad
#'   \bm{X}_j = \bigl[\bm{1},\; \hat{\bm{\tau}}_{A_j}\bigr]
#' }
#'
#' where \eqn{\hat{\bm{\tau}}_j} and \eqn{\hat{\bm{\tau}}_{A_j}} are the
#' \eqn{n}-vectors of predicted values (BLUEs) for treatment \eqn{j} and its
#' conditioning set across genotypes.  The **response index**
#' \eqn{\tilde{\bm{\tau}}_j = \bar{\bm{H}}_j \hat{\bm{\tau}}_j} is the
#' vector of OLS residuals, where
#' \eqn{\bar{\bm{H}}_j = \bm{I} - \bm{X}_j(\bm{X}_j^\top\bm{X}_j)^{-1}\bm{X}_j^\top}
#' is the hat-matrix annihilator.  Its variance--covariance matrix is
#' \eqn{\sigma_j^2\bar{\bm{H}}_j} with \eqn{df = n - a - 1} degrees of
#' freedom.
#'
#' The same four conditioning schemes as [randomRegress()] are available
#' via the `type` argument:
#'
#' \describe{
#'   \item{`"baseline"` (default)}{Each non-first treatment regressed on
#'     \code{levs[1]} alone (simple regression, \eqn{df = n-2}).}
#'   \item{`"sequential"`}{Treatment \eqn{j} regressed on all preceding
#'     treatments \code{levs[1:(j-1)]} (multiple regression).  By the
#'     Gram--Schmidt property, residuals from treatment \eqn{j} are
#'     uncorrelated with all previous treatment BLUEs.}
#'   \item{`"partial"`}{Each treatment regressed on all other treatments
#'     simultaneously.  Residuals are partial regression residuals,
#'     uncorrelated with the conditioning set but not necessarily with each
#'     other.}
#'   \item{`"custom"`}{Conditioning set for each treatment supplied
#'     explicitly via `cond`.}
#' }
#'
#' @param model   An ASReml-R V4 model object.
#' @param term    Character string specifying the `classify` term containing
#'   the treatment and genotype factors, e.g. `"Treatment:Genotype"` or
#'   `"Treatment:Site:Genotype"`.
#' @param by      Optional character string naming a factor in `term` to split
#'   the analysis by (e.g. `"Site"`).  Separate regressions are performed
#'   within each level of `by`.  When `NULL` (default) the regression is
#'   performed over all observations jointly.
#' @param levs    Character vector of length \eqn{\ge 2} giving the treatment
#'   labels.  The first element is the baseline for `type = "baseline"` and
#'   `type = "sequential"`.  Cannot be `NULL`.
#' @param type    Conditioning scheme.  One of `"baseline"` (default),
#'   `"sequential"`, `"partial"`, or `"custom"`.  See **Description**.
#' @param cond    Named list required when `type = "custom"`.  Each name must
#'   be a treatment label from `levs`; each value is `NULL` (unconditional)
#'   or a character vector of conditioning treatment labels.  Treatments
#'   absent from `cond` are treated as unconditional.
#' @param min_obs Minimum number of genotypes with non-missing BLUEs in all
#'   required treatments for a regression to be attempted.  Defaults to
#'   \eqn{\max(5,\; 2(|A_j|_{\max} + 1))}.  Groups below this threshold
#'   produce a warning and are omitted.
#' @param ...     Additional arguments forwarded to
#'   `asreml::predict.asreml()`.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{`blues`}{Data frame with columns: the grouping variable (or
#'     `"Group"` when `by = NULL`), `Genotype`, one raw BLUE column per
#'     element of `levs`, one `resp.<lev>` column per conditioned treatment
#'     (response indices / OLS residuals), one `se.<lev>` column (residual
#'     standard errors), and one `HSD.<lev>` column (Tukey HSD for pairwise
#'     genotype comparisons on the response-index scale).}
#'   \item{`beta`}{Named list of length \eqn{n_{\text{cond}}}.  Each element
#'     `beta[["<lev>"]]` is a data frame with the grouping variable plus one
#'     column per conditioning treatment, giving the OLS regression
#'     coefficients estimated within each group.}
#'   \item{`sigmat`}{Numeric matrix of dimensions
#'     \eqn{n_{\text{groups}} \times n_{\text{cond}}} containing the residual
#'     standard deviation \eqn{\sigma_j} from each OLS regression.}
#'   \item{`cond_list`}{The resolved conditioning structure (named list).}
#'   \item{`type`}{The `type` argument used.}
#' }
#'
#' @seealso [randomRegress()], `asreml::predict.asreml()`
#'
#' @examples
#' \dontrun{
#' ## Baseline: T1 and T2 each regressed on T0 within each Site
#' res_base <- fixedRegress(model,
#'                             term = "Treatment:Site:Genotype",
#'                             by   = "Site",
#'                             levs = c("T0","T1","T2"))
#'
#' ## Sequential: T1|T0, T2|{T0,T1} -- Gram-Schmidt orthogonal residuals
#' res_seq  <- fixedRegress(model,
#'                             term = "Treatment:Site:Genotype",
#'                             by   = "Site",
#'                             levs = c("T0","T1","T2"),
#'                             type = "sequential")
#'
#' ## Partial: each treatment vs all others within each site
#' res_part <- fixedRegress(model,
#'                             term = "Treatment:Site:Genotype",
#'                             by   = "Site",
#'                             levs = c("T0","T1","T2"),
#'                             type = "partial")
#'
#' ## Custom: T0 unconditional, T1|T0, T2|{T0,T1}
#' res_cust <- fixedRegress(model,
#'                             term = "Treatment:Site:Genotype",
#'                             by   = "Site",
#'                             levs = c("T0","T1","T2"),
#'                             type = "custom",
#'                             cond = list(T0 = NULL,
#'                                         T1 = "T0",
#'                                         T2 = c("T0","T1")))
#' }
#'
#' @export
fixedRegress <- function(model, term = "Treatment:Genotype",
                            by   = NULL, levs = NULL,
                            type = "baseline", cond = NULL,
                            min_obs = NULL, ...) {

  # ---- Validate inputs and build conditioning structure ------------------
  if (is.null(levs) || length(levs) < 2L)
    stop("At least two treatment levels must be supplied in 'levs'.")

  ntreat <- length(levs)
  type   <- match.arg(type, c("baseline", "sequential", "partial", "custom"))

  cond_list   <- .condList(levs, type, cond)
  conditioned <- levs[!vapply(cond_list, is.null, logical(1L))]
  n_cond      <- length(conditioned)
  if (n_cond == 0L)
    stop("No treatments have a conditioning set. Check 'type' or 'cond'.")

  # ---- Parse term --------------------------------------------------------
  pterm <- term
  terms <- unlist(strsplit(term, ":"))
  if (length(terms) < 2L)
    stop("'term' must contain at least two ':'-separated variables.")

  # ---- Get BLUEs ---------------------------------------------------------
  pred <- predict(model, classify = pterm, ...)
  whna <- !is.na(pred$pvals$predicted.value)
  pv   <- pred$pvals[whna, ]

  # Identify the treatment variable (factor whose levels contain all of levs)
  has_levs <- vapply(terms, function(v) {
    col <- pv[[v]]
    is.factor(col) && all(levs %in% levels(col))
  }, logical(1L))

  if (!any(has_levs))
    stop("Treatment levels in 'levs' not found in any variable of 'term'.")
  tnam <- terms[has_levs][1L]

  # ---- Identify the grouping variable (by) and regression variable -------
  if (!is.null(by)) {
    bys <- unlist(strsplit(by, ":"))
    if (!all(bys %in% terms))
      stop("Variables in 'by' are not all present in 'term'.")
    if (tnam %in% bys)
      stop("Treatment variable '", tnam, "' cannot also be in 'by'.")
    rterm <- terms[!(terms %in% c(tnam, bys))]
    if (!length(rterm))
      stop("No regression variable remains after removing treatment and 'by'.")
    if (length(bys) > 1L)
      pv[["by__"]] <- apply(pv[, bys, drop = FALSE], 1L,
                             paste, collapse = ":")
    else
      pv[["by__"]] <- as.character(pv[[bys[1L]]])
    by_col <- by
  } else {
    rterm  <- terms[terms != tnam]
    pv[["by__"]] <- "All"
    by_col <- "Group"
  }

  # Create regression key (identifies genotype/unit within a group)
  if (length(rterm) > 1L)
    pv[["reg__"]] <- apply(pv[, rterm, drop = FALSE], 1L, paste, collapse = ":")
  else
    pv[["reg__"]] <- as.character(pv[[rterm[1L]]])

  gnam  <- paste(rterm, collapse = ":")   # genotype column name in output
  um    <- unique(pv[["by__"]])
  ns    <- length(um)

  # ---- Default minimum observations --------------------------------------
  max_a <- max(vapply(cond_list[conditioned], length, integer(1L)))
  if (is.null(min_obs))
    min_obs <- max(5L, 2L * (max_a + 1L))

  # ---- Column name vectors -----------------------------------------------
  resp_nams <- paste0("resp.", conditioned)
  se_nams   <- paste0("se.",   conditioned)
  hsd_nams  <- paste0("HSD.",  conditioned)

  # ---- Initialise beta storage -------------------------------------------
  # beta[[lv_j]]: data frame with by_col + one column per conditioning treatment
  beta_list <- setNames(vector("list", n_cond), conditioned)
  for (lv in conditioned) {
    A_j             <- cond_list[[lv]]
    beta_list[[lv]] <- as.data.frame(
      matrix(NA_real_, ns, length(A_j),
             dimnames = list(NULL, A_j))
    )
    beta_list[[lv]] <- cbind(
      setNames(data.frame(um, stringsAsFactors = FALSE), by_col),
      beta_list[[lv]]
    )
  }

  sigmat <- matrix(NA_real_, ns, n_cond, dimnames = list(um, conditioned))

  # ---- Main loop: one iteration per group --------------------------------
  result_list <- vector("list", ns)

  for (i in seq_along(um)) {

    pvg <- pv[pv[["by__"]] == um[i], ]

    # Wide matrix of BLUEs: rows = genotypes, cols = treatments
    # Only keep genotypes present in every required treatment
    geno_by_trt <- lapply(levs, function(lv) {
      rows <- pvg[pvg[[tnam]] == lv, ]
      setNames(rows$predicted.value, rows[["reg__"]])
    })
    names(geno_by_trt) <- levs

    # Common genotypes across ALL treatment levels
    common_genos <- Reduce(intersect, lapply(geno_by_trt, names))

    if (length(common_genos) < min_obs) {
      warning("Group '", um[i], "' has only ", length(common_genos),
              " common genotypes (< min_obs = ", min_obs, "). Skipping.")
      next
    }

    n <- length(common_genos)

    # Aligned BLUE matrix [n x ntreat]
    blue_mat <- matrix(NA_real_, n, ntreat,
                       dimnames = list(common_genos, levs))
    for (lv in levs)
      blue_mat[, lv] <- geno_by_trt[[lv]][common_genos]

    # Output matrices for this group [n x n_cond]
    resp_mat <- matrix(NA_real_, n, n_cond, dimnames = list(common_genos, resp_nams))
    se_mat   <- matrix(NA_real_, n, n_cond, dimnames = list(common_genos, se_nams))
    hsd_mat  <- matrix(NA_real_, n, n_cond, dimnames = list(common_genos, hsd_nams))

    # ---- Per-treatment regression ----------------------------------------
    for (ci in seq_len(n_cond)) {

      lv_j <- conditioned[ci]
      A_j  <- cond_list[[lv_j]]
      a    <- length(A_j)
      ndf  <- n - a - 1L

      if (ndf < 3L) {
        warning("Group '", um[i], "' treatment '", lv_j, "': df = ", ndf,
                " (< 3). Skipping.")
        next
      }

      blues_j  <- blue_mat[, lv_j]
      blues_Aj <- blue_mat[, A_j, drop = FALSE]

      # OLS regression: blues_j ~ blues_Aj
      fit <- tryCatch(
        lm(blues_j ~ blues_Aj),
        error = function(e) {
          warning("OLS failed for group '", um[i], "' treatment '",
                  lv_j, "': ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit)) next

      resp <- residuals(fit)
      sig  <- summary(fit)$sigma

      # Hat-matrix annihilator: H-bar = I - X(X'X)^{-1}X'
      Xm   <- model.matrix(fit)
      XtX_inv <- tryCatch(
        solve(crossprod(Xm)),
        error = function(e) {
          warning("Near-singular design matrix for group '", um[i],
                  "' treatment '", lv_j, "'. Skipping.")
          NULL
        }
      )
      if (is.null(XtX_inv)) next

      Hbar <- diag(n) - Xm %*% tcrossprod(XtX_inv, Xm)
      rv   <- sig^2 * Hbar                           # residual var-cov

      # SED matrix and Tukey HSD
      dv   <- diag(rv)
      sed  <- outer(dv, dv, "+") - 2 * rv
      sed  <- sed[lower.tri(sed)]
      sed[sed < 0] <- NA_real_
      hsd  <- (mean(sqrt(sed), na.rm = TRUE) / sqrt(2)) *
                qtukey(0.95, n, df = ndf)

      resp_mat[, ci] <- resp
      se_mat[, ci]   <- sqrt(dv)
      hsd_mat[, ci]  <- hsd
      sigmat[i, ci]  <- sig

      # Store OLS coefficients (drop intercept, name by conditioning treatment)
      bcoef         <- coef(fit)[-1L]
      names(bcoef)  <- A_j
      beta_list[[lv_j]][i, A_j] <- bcoef
    }

    # ---- Assemble group result -------------------------------------------
    result_list[[i]] <- cbind(
      setNames(data.frame(rep(um[i], n), common_genos,
                          stringsAsFactors = FALSE),
               c(by_col, gnam)),
      as.data.frame(blue_mat),
      as.data.frame(resp_mat),
      as.data.frame(se_mat),
      as.data.frame(hsd_mat)
    )
  }

  # ---- Assemble final output ---------------------------------------------
  result_list <- result_list[!vapply(result_list, is.null, logical(1L))]
  blues_out   <- do.call(rbind, result_list)
  rownames(blues_out) <- NULL

  list(blues     = blues_out,
       beta      = beta_list,
       sigmat    = sigmat,
       cond_list = cond_list,
       type      = type)
}
